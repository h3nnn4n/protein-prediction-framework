from mpi4py import MPI
import config_loader
import results
import random
import comm
import sys
import os
import de


if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    size = comm.size
    rank = comm.rank
    status = MPI.Status()

    if rank == 0:
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
        random.shuffle(todo)

        ##################

        sent = 0
        recived = 0
        alive = list(range(1, size))

        for i in alive:
            comm.send(('wakeup', None), dest=i)

        while recived < len(todo):
            s, d = comm.recv(source=MPI.ANY_SOURCE, status=status)

            if s == 'ok':
                recived += 1
                # print(s, d, sent, recived, len(todo), status.Get_source())
                d.store()

            if sent < len(todo):
                cf = config_loader.ConfigLoader(todo[sent])
                comm.send(('ok', cf), dest=status.Get_source())
                sent += 1
            else:
                print('sent die to %d' % (status.Get_source()))
                comm.send(('die', None), dest=status.Get_source())
                alive.remove(status.Get_source())

        for i in alive:
            comm.send(('die', None), dest=status.Get_source())

        print('master finished')
    else:
        comm.send(('boot', 0), dest=0)
        s = ''
        jobs_finished = 0
        while s != 'die':
            s, cf = comm.recv(source=MPI.ANY_SOURCE, status=status)

            if s == 'ok':
                jobs_finished += 1

                c = None  # comm.Pigeon()
                d = None

                d = de.DE(pop_size=cf['pop_size'], max_iters=cf['max_iters'], pname=cf['pname'],
                          f_factor=cf['f_factor'], c_rate=cf['c_rate'], allatom=cf['allatom'])

                d.set_coms(c)

                cf.inject(d)

                d.run()

                comm.send(('ok', results.Results(d.stats_name, d.config_name)), dest=status.Get_source())

        print(' %2d finished %3d jobs' % (rank, jobs_finished))
