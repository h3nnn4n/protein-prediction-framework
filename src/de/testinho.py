from rosetta_pack import RosettaPack
import pyrosetta
import random


rp = RosettaPack(name='1zdd')

temp = 2.0
n_moves = 5

allatom = False
pose = pyrosetta.pose_from_sequence(rp.target)

if allatom:
    score_function = pyrosetta.get_fa_scorefxn()
else:
    centroid_switch = pyrosetta.SwitchResidueTypeSetMover('centroid')
    score_function = pyrosetta.create_score_function('score3')
    centroid_switch.apply(pose)

mc = pyrosetta.MonteCarlo(pose, score_function, temp)

fragset3 = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(3)
fragset9 = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(9)

fragset3.read_fragment_file('../../protein_data/1zdd/output/1zdd.200.3mers')
fragset9.read_fragment_file('../../protein_data/1zdd/output/1zdd.200.9mers')

cost = pyrosetta.rosetta.protocols.simple_moves.GunnCost()

movemap = pyrosetta.MoveMap()
movemap.set_bb(True)
movemap = rp.get_new_movemap_with_free_coil_and_loop()

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

if not allatom:
    minmover.score_function(score_function)

perturbation_seq = pyrosetta.rosetta.protocols.moves.SequenceMover()
perturbation_seq.add_mover(smallmover)
perturbation_seq.add_mover(shearmover)
perturbation_trial = pyrosetta.TrialMover(perturbation_seq, mc)
perturbation_rep = pyrosetta.rosetta.protocols.moves.RepeatMover(perturbation_trial, 10)

frag_trial_9mer = pyrosetta.TrialMover(mover_9mer, mc)
frag_trial_3mer_smooth = pyrosetta.TrialMover(mover_3mer_smooth, mc)
frag_sequence = pyrosetta.rosetta.protocols.moves.SequenceMover()
frag_sequence.add_mover(frag_trial_9mer)
frag_sequence.add_mover(frag_trial_3mer_smooth)
frag_repeater = pyrosetta.rosetta.protocols.moves.RepeatMover(frag_sequence, 100)

hopper = pyrosetta.rosetta.protocols.moves.SequenceMover()
hopper.add_mover(perturbation_rep)
hopper.add_mover(frag_repeater)
hopper.add_mover(minmover)

hopper_trial = pyrosetta.TrialMover(hopper, mc)
hopper_rep = pyrosetta.rosetta.protocols.moves.RepeatMover(hopper_trial, 5)


def set_pose(pose, data):
    for i in range(len(rp.target)):
        pose.set_phi(i + 1, data[i * 3 + 0])
        pose.set_psi(i + 1, data[i * 3 + 1])
        pose.set_omega(i + 1, data[i * 3 + 2])


def reset_pose():
    angles = [random.uniform(-180, 180) for _ in range(3 * len(rp.target))]
    set_pose(pose, angles)
    score = score_function(pose)
    mc.set_last_accepted(score)
    mc.set_last_accepted_pose(pose)
    mc.reset_counters()


def run(evals=500000):
    print("%8d %8d %12.4f %12.4f" % (0, mc.total_trials(), score_function(pose), mc.lowest_score()))
    i = 0
    while True:
        hopper_rep.apply(pose)
        print("%8d %8d %12.4f %12.4f" % (i + 1, mc.total_trials(), score_function(pose), mc.lowest_score()))
        if random.random() < 0.05:
            mc.recover_low(pose)

        if mc.total_trials() >= evals:
            break

        i += 1
