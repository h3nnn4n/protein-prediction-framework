from mpi4py import MPI


class Pigeon:
    def __init__(self):
        self.comm = MPI.COMM_WORLD
        self.size = self.comm.size
        self.rank = self.comm.rank
        self.status = MPI.Status()

        print('%3d is on %s' % (self.rank, MPI.Get_processor_name()))

    def migration(self, best):
        if self.size <= 1:
            return None

        debug = False
        # debug = True

        comm = self.comm

        if self.rank == 0:
            dest = 1
            source = self.size - 1

            if debug:
                print('%d sending to %d' % (self.rank, dest))

            comm.send(best, dest=dest)
            r = comm.recv(source=source, status=self.status)

            if debug:
                print("%d got message from %d" % (self.rank, self.status.Get_source()))

            return r
        else:
            source = (self.rank - 1)
            dest = (self.rank + 1) % (self.size)

            if debug:
                print('%d sending to %d' % (self.rank, dest))

            r = comm.recv(source=source)
            comm.send(best, dest=dest)

            if debug:
                print("%d got message from %d" % (self.rank, self.status.Get_source()))

            return r
