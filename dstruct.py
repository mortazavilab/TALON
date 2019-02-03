class Struct(dict):
    """
    Make a dict behave as a struct.

    Example:
    
        test = Struct(a=1, b=2, c=3)

    """
    def __init__(self,**kw):
        dict.__init__(self,kw)
        self.__dict__ = self
