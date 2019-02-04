import multiprocessing
import sys
import os
from rosetta_classic_abinitio import ClassicAbinitio


def wrap(p):
        pname, factor = p
        ClassicAbinitio(pname).run(factor)


def main():
    np = 4

    conf, reps = get_conf_and_reps()
    todo = open_todo_and_get_todolist(conf)

    todo *= reps

    p = multiprocessing.Pool(np)
    p.map(wrap, todo)


def get_conf_and_reps():
    reps = 1
    if len(sys.argv) > 1:
        for arg in sys.argv[1:]:
            if os.path.isfile(arg):
                conf = arg
            else:
                reps = int(arg)
    else:
        raise NotImplementedError('Running with not args is not supported')

    return conf, reps


def open_todo_and_get_todolist(conf):
    todo = []
    with open(conf) as f:
        for l in f.readlines():
            pname, factor = l.strip().split(' ')
            todo.append((pname, int(factor)))

    return todo


if __name__ == '__main__':
    main()
