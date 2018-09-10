import pyrosetta


target = 'MTKQEKTALNMARFIRSQTLTLLEKLNELDADEQADICESLHDHADELYRSCLARFGDDGENL'
temp = 1.0
n_moves = 1


pyrosetta.init('-out:level 0')

scorefxn = pyrosetta.get_fa_scorefxn()
pose = pyrosetta.pose_from_sequence(target)

movemap = pyrosetta.MoveMap()
movemap.set_bb(True)

smallmover = pyrosetta.rosetta.protocols.simple_moves.SmallMover(movemap, temp, n_moves)
shearmover = pyrosetta.rosetta.protocols.simple_moves.ShearMover(movemap, temp, n_moves)
minmover = pyrosetta.rosetta.protocols.simple_moves.MinMover()

mc = pyrosetta.MonteCarlo(pose, scorefxn, 2.0)
trial = pyrosetta.TrialMover(smallmover, mc)

for i in range(100):
    trial.apply(pose)
    print(pose.energies().total_energy())
