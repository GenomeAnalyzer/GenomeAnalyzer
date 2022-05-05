def read_file(name):
    f = open(name, "r")
    fr = f.readlines()
    return "".join([f[:-1] for f in fr[1:]])
