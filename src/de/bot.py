import multiprocessing
import main
import sys
import os


if __name__ == '__main__':
    conf = None
    pwd = None
    reps = 1
    todo = []
    np = 8

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

    # print('running with %s %s %s' % (conf, pwd, reps))

    with open(conf) as f:
        for l in f.readlines():
            line = l.rstrip().lstrip()
            if os.path.isfile(line):
                todo.append(line)
            else:
                print('%s was not found, skiping.' % (line))

    todo *= reps

    def wrap(p):
        # print('starting at %s' % (os.getcwd()))
        os.chdir(pwd)
        main.boot(p)
        # print('finishing at %s' % (os.getcwd()))

    p = multiprocessing.Pool(np)
    p.map(wrap, todo)
