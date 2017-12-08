import config_loader
import protein_data
import rosetta_pack
import comm
import sys
import de


if __name__ == '__main__':
    c = comm.Pigeon()
    d = None

    if len(sys.argv) > 1:
        conf_file = sys.argv[1]
    else:
        conf_file = None

    cf = config_loader.ConfigLoader(conf_file)

    if c.rank == 0:
        d = de.DE(pop_size=cf['pop_size'], max_iters=cf['max_iters'], pname=cf['pname'],
                f_factor=cf['f_factor'], c_rate=cf['c_rate'], allatom=cf['allatom'])

        d.set_coms(c)

        cf.inject(d)

        d.run()
    else:
        rosetta_pack_ = rosetta_pack.RosettaPack(cf['pname'])
        pdata = protein_data.ProteinData(rosetta_pack_, allatom=cf['allatom'])

        for i in range(cf['max_iters']):
            # print('new gen on %d' % c.rank)
            c.comm.send((-1, 0), dest=0)
            code = 0
            while code >= 0:
                code, angles = c.comm.recv(source=0)

                if code >= 0:
                    pdata.new_angles(angles)
                    pdata.eval()

                    c.comm.send((code, pdata.score), dest=0)
            # print('%d barrier' % c.rank)
            c.comm.barrier()
