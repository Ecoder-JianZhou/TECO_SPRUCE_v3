

def a():
    aa = 1
    def b():
        bb = 2
        res = aa + bb
        return res
    return b
print(a()()) 