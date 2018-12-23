import pyrosetta


target = 'MTKQEKTALNMARFIRSQTLTLLEKLNELDADEQADICESLHDHADELYRSCLARFGDDGENL'
temp = 1.0
n_moves = 1


pyrosetta.init('-out:level 0')

scorefxn = pyrosetta.get_fa_scorefxn()
pose = pyrosetta.pose_from_sequence(target)
mc = pyrosetta.MonteCarlo(pose, scorefxn, 2.0)

movemap = pyrosetta.MoveMap()
movemap.set_bb(True)

smallmover = pyrosetta.rosetta.protocols.simple_moves.SmallMover(movemap, temp, n_moves)
shearmover = pyrosetta.rosetta.protocols.simple_moves.ShearMover(movemap, temp, n_moves)
minmover = pyrosetta.rosetta.protocols.simple_moves.MinMover()

perturbation_seq = pyrosetta.rosetta.protocols.moves.SequenceMover()
perturbation_seq.add_mover(smallmover)
perturbation_seq.add_mover(shearmover)
perturbation_rep = pyrosetta.rosetta.protocols.moves.RepeatMover(perturbation_seq, 10)

hopper = pyrosetta.rosetta.protocols.moves.SequenceMover()
hopper.add_mover(perturbation_rep)
hopper.add_mover(minmover)

hopper_trial = pyrosetta.TrialMover(hopper, mc)
hopper_rep = pyrosetta.rosetta.protocols.moves.RepeatMover(perturbation_seq, 100)

hopper_rep.apply(pose)
