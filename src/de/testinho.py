import pyrosetta


target = 'MTKQEKTALNMARFIRSQTLTLLEKLNELDADEQADICESLHDHADELYRSCLARFGDDGENL'
temp = 1.0
n_moves = 5

pyrosetta.init('-out:level 0')

scorefxn = pyrosetta.get_fa_scorefxn()
pose = pyrosetta.pose_from_sequence(target)
mc = pyrosetta.MonteCarlo(pose, scorefxn, temp)

movemap = pyrosetta.MoveMap()
movemap.set_bb(True)

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
