import numpy as np
from enum import Enum


def checktype(x, xtype):
    """Perform type checking.
    Params:
    ______
    x : atomic, iterable.
        If atomic, xtype must be a type, otherwise xtype is assusmed to be an iterbable of
        the same size as x.
    xtype : atomic, iterable, iterable of tuples
        If atomic, x must be an `instance`, otherwise x is assumed to be list of instances.
        
    example :
    ________
    checktype(1,int)
    checktype([1, 2], [int, (int, str)])
        """
    if not isinstance(xtype, type):
        for val, valtype in zip(x, xtype) :
            if not isinstance(valtype, type) :
                at_least_one_type = False
                for ith_type in valtype :
                    assert isinstance(ith_type, type)
                    at_least_one_type = isinstance(val, ith_type)
                    if at_least_one_type :
                        break
                assert at_least_one_type
            else :
                assert isinstance(val, valtype)
    else :
        assert isinstance(x, xtype)
        
class btype(Enum):
    UTILITY = "U"
    U = "U" # An alias for UTILITY
    RESIDENCE = "R"
    R = "R" # An alias for RESIDENCE
    def is_utility(self):
        return self.value == "U"
    def is_residence(self):
        return self.value == "R"      
        
        
        
#%%
def project_from_file(file):
    with open(file, "r") as f :
        HWDB = f.readline()
        H,W,D,B = map(int, HWDB.replace("\n", "").split())
        
        bdgs = []
        for b in range(B):
            bdg = {}
            thwur = f.readline()
            h,w,ur = map(int, thwur[1:].replace("\n", "").split())
            t = thwur[0]
            if t == "R":
                bdg["kind"] = btype.RESIDENCE
                bdg["capacity"] = ur
                bdg["utype"] = None
                
            elif t == "U":
                bdg["kind"] = btype.UTILITY
                bdg["capacity"] = None
                bdg["utype"] = ur

            else :
                raise ValueError("A value different from 'R' and 'U' is found as project value")
            
            grid = []
            for i in range(h):
                row = f.readline()
                row = [ 1 if c == "#" else 0 for 
                        c in row.replace("\n", "") ]
                grid.append(row)
            # print("\n", grid)
            bdg["grid"] = np.array(grid, dtype = int)
            bdgs.append(bdg)
        return {"city_params": {"H":H, "W":W, "D":D, "B":B},"city_buildings" : bdgs} 