
#%%
import numpy as np
from enum import Enum
import itertools


from utils import project_from_file, checktype, btype

from hyperopt import fmin, tpe, hp,  STATUS_OK, Trials
from matplotlib import pyplot as plt


#%%
class plan(object):
    def __init__(self, H=2, W=2, grid = None):
        
        checktype([H, W], [int, int])
        
        if grid is None :
            grid = np.zeros((H,W))
        else :
            if not isinstance(grid, np.ndarray) :
                grid = np.array(grid, dtype = int)
            H,W = grid.shape
        self.grid = grid
        self.H = H
        self.W = W
    
    def reset(self):
        self.grid = np.zeros(self.grid.shape)
    
    def merge(self, top_left_edge_pos, other):
        i0 = top_left_edge_pos[0]
        j0 = top_left_edge_pos[1]
        i1 = i0 + other.grid.shape[0]
        j1 = j0 + other.grid.shape[1]
        v = self.grid[i0:i1, j0:j1 ][other.grid != 0]
        feedback = {"built": False, "i0":i0, "i1":i1, "j0":j0, "j1":j1}
        if not v.any():
            self.grid[i0:i1, j0:j1 ] += other.grid
            feedback["built"] = True
        return feedback

        
#%%
def manhattan(cell1, cell2):
    return np.abs(cell2[0] - cell1[0]) + np.abs(cell2[1] - cell1[1])


#%%
class building(plan):
    def __init__(self,H=2, W=2, grid=None, kind = btype.RESIDENCE, capacity= None, utype = None, name = "BuildingX"):
        checktype([kind, capacity, utype, name], 
                  [ btype, (type(None),int), (type(None),int), str])
        super().__init__(H,W, grid)
        self.kind = kind
        self.capacity = capacity
        self.utype = utype
        self.name = name
        
    def is_residence(self):
        return self.kind.is_residence()
    
    def is_utility(self):
        return self.kind.is_utility()



#%%
class city(plan):
    def __init__(self, H=2,W=2,grid=None, D = 5,B = 0, name = "cityX"):
        checktype([name], [ str])
        super().__init__(H,W, grid)
        self.name = name
        self.D = D
        self.B = B
        self.residences = []
        self.utilities = []
        self.manhattans = []
        self.score = None
        self.score_details = []
        self.to_builds = []
        
    def reset(self):
        self.grid = np.zeros(self.grid.shape)
        self.residences = []
        self.utilities = []
        self.manhattans = []
        self.score = None
        self.score_details = []
        
    def build(self, top_left_pos, to_build):
        feedback = self.merge(top_left_pos, to_build)
        if feedback["built"]:
            feedback["name"] = to_build.name
            
#             feedback["kind"] = to_build.kind
            if  to_build.is_residence() :
                feedback["capacity"] = to_build.capacity
                self.residences.append(feedback)
            elif  to_build.is_utility() :
                feedback["utype"] = to_build.utype
                self.utilities.append(feedback)
            else :
                raise ValueError("The building kind must be %s "% btype)
            
#             self.B += 1
            
        return feedback["built"]
            
            
    def get_building_cells(self, index, kind = btype.RESIDENCE):
        if kind.is_residence() :
            b = self.residences[index]
        else : 
            b = self.utilities[index]
        i0,i1,j0,j1 = b["i0"], b["i1"], b["j0"], b["j1"]
        return itertools.product(np.arange(i0,i1), np.arange( j0, j1 ) )
    
    def distance(self, rindex, uindex):
        for dist in self.manhattans :
            keys= dist.keys()
            if rindex in keys and uindex in keys :
                return dist
        
        minmanh = +np.inf
        bcells = self.get_building_cells(rindex, kind = btype.RESIDENCE)
        ucells = self.get_building_cells(uindex, kind = btype.UTILITY)
        
        for bcell, ucell in itertools.product(bcells,ucells) :
            manh = manhattan(bcell, ucell)
            if manh < minmanh:
                minmanh = manh
                minbcell = bcell
                minucell = ucell
            if minmanh == 0 :
                break
        self.manhattans.append({"rindex": rindex, "uindex": uindex, "bcell": minbcell, "ucell": minucell,
                                "manh": minmanh})
        return self.manhattans[-1]
    
    def scorer(self):
        s = 0
        for r in range(len(self.residences)) : 
            utility_types = set()
            for u in range(len(self.utilities)) :
            
                if self.distance(r, u)["manh"] <= self.D :
                    capacity = self.residences[r]["capacity"]
                    utype = self.utilities[u]["utype"]
                    gain = capacity*(not utype in utility_types)
                    s += gain
                    utility_types.add(utype)
                    self.score_details.append({"rindex": r, "uindex": u, "gain": gain, "utype": utype })
        self.score = s
        return self.score
    def set_tobuilds(self,to_builds):
        self.to_builds = to_builds
        self.B = len(to_builds)
        
#     def set_tobuilds_from_file(self, file):
#         self.set_tobuilds(self.city_from_file(file))
    
    def builder__scorer(self, top_left_pos, to_builds = None):
        self.reset()
        if to_builds is None :
            to_builds = self.to_builds
        for pos, to_build in zip(top_left_pos, to_builds) :
#             assert all([pos[0]>=0, pos[1] >= 0, pos[0] < self.grid.H, pos[1] < self.grid.W])
            try :
                pos = (abs(int(pos[0])), int(abs(pos[1])) )
                assert all([ pos[0] < self.H, pos[1] < self.W])
                self.build(pos, to_build)
            except :
                pass
        try:
            score = -self.scorer()
        except :
            score = +np.Inf
        
        return {"loss": score, "status":  STATUS_OK }
    
    def project_from_file(file):
        return project_from_file(file)
        
    
    def city_from_file(file):
        pjt = project_from_file(file)
        H, W, D = (pjt["city_params"]["H"],pjt["city_params"]["W"]
                    ,pjt["city_params"]["D"])
        c = city(H,W, D= D)

        bdgs = pjt["city_buildings"]
        for i in range(len(bdgs)) :
            bdg = bdgs[i]
            bdg = building(0,0,bdg["grid"], bdg["kind"], bdg["capacity"], bdg["utype"], "Building%d"%i)

            bdgs[i] = bdg
        c.set_tobuilds(bdgs)
        return c
