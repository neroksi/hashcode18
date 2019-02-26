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