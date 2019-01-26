import pyrosetta


target = 'MTKQEKTALNMARFIRSQTLTLLEKLNELDADEQADICESLHDHADELYRSCLARFGDDGENL'
temp = 1.0
n_moves = 5

pyrosetta.init('-out:level 0')

scorefxn = pyrosetta.get_fa_scorefxn()
pose = pyrosetta.pose_from_sequence(target)
mc = pyrosetta.MonteCarlo(pose, scorefxn, temp)

fragset3 = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(3)
fragset9 = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(9)

fragset3.read_fragment_file('../../protein_data/1zdd/output/1zdd.200.3mers')
fragset9.read_fragment_file('../../protein_data/1zdd/output/1zdd.200.9mers')

cost = pyrosetta.rosetta.protocols.simple_moves.GunnCost()

movemap = pyrosetta.MoveMap()
movemap.set_bb(True)

mover_3mer = pyrosetta.rosetta.protocols.simple_moves.ClassicFragmentMover(fragset3, movemap)
mover_9mer = pyrosetta.rosetta.protocols.simple_moves.ClassicFragmentMover(fragset9, movemap)
mover_3mer_smooth = pyrosetta.rosetta.protocols.simple_moves.SmoothFragmentMover(fragset3, movemap, cost)
mover_9mer_smooth = pyrosetta.rosetta.protocols.simple_moves.SmoothFragmentMover(fragset9, movemap, cost)

mover_3mer.set_movemap(movemap)
mover_9mer.set_movemap(movemap)
mover_3mer_smooth.set_movemap(movemap)
mover_9mer_smooth.set_movemap(movemap)

smallmover = pyrosetta.rosetta.protocols.simple_moves.SmallMover(movemap, temp, n_moves)
shearmover = pyrosetta.rosetta.protocols.simple_moves.ShearMover(movemap, temp, n_moves)
minmover = pyrosetta.rosetta.protocols.minimization_packing.MinMover()

perturbation_seq = pyrosetta.rosetta.protocols.moves.SequenceMover()
perturbation_seq.add_mover(smallmover)
perturbation_seq.add_mover(shearmover)
perturbation_rep = pyrosetta.rosetta.protocols.moves.RepeatMover(perturbation_seq, 10)

hopper = pyrosetta.rosetta.protocols.moves.SequenceMover()
hopper.add_mover(perturbation_rep)
hopper.add_mover(minmover)

hopper_trial = pyrosetta.TrialMover(hopper, mc)
hopper_rep = pyrosetta.rosetta.protocols.moves.RepeatMover(perturbation_seq, 100)


def run(n=100):
    for i in range(n):
        hopper_rep.apply(pose)
        print("%8d %12.4f" % (i, scorefxn(pose)))
