import multiprocessing
import sys
import os
from main import boot


def wrap(p):
    print(p)
    boot(p)


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
            line = l.strip()
            if os.path.isfile(line):
                todo.append(line)
            else:
                print('%s was not found, skiping.' % (line))

    return todo


if __name__ == '__main__':
    main()
