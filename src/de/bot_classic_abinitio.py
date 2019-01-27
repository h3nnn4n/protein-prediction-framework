import multiprocessing
import sys
import os

from rosetta_classic_abinitio import ClassicAbinitio

if __name__ == '__main__':
    conf = None
    pwd = None
    reps = 1
    todo = []
    np = 4

    if len(sys.argv) > 1:
        for arg in sys.argv[1:]:
            if os.path.isfile(arg):
                conf = arg
            elif os.path.isdir(arg):
                pwd = arg
            else:
                reps = int(arg)

        if pwd is None:
            pwd = os.getcwd()
    else:
        sys.exit()

    print('running with %s %s %s' % (conf, pwd, reps))

    with open(conf) as f:
        for l in f.readlines():
            pname, factor = l.strip().split(' ')
            todo.append((pname, int(factor)))

    todo *= reps

    def wrap(p):
        # print('running with params: %s %s' % p)
        # print('starting at %s' % (os.getcwd()))
        pname, factor = p
        ClassicAbinitio(pname).run(factor)
        # print('finishing at %s' % (os.getcwd()))

    p = multiprocessing.Pool(np)
    p.map(wrap, todo)
