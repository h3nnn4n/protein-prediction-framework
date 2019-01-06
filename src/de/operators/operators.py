import random
import types
from operators.rand1bin_global import rand1bin_global


class Operators:
    def __init__(self, de=None):
        self.de = de
        self.available_operators = {
            'rand1bin_global': self.get_rand1bin_global
        }

    def get_random_individual(self):
        p = random.randint(0, self.de.pop_size - 1)
        return self.de.pop[p]

    def get_random_individuals(self, n=1):
        return random.sample(self.de.pop, n)

    def get_operator(self, name):
        return self.available_operators[name]()

    def get_rand1bin_global(self):
        self.rand1bin_global = types.MethodType(rand1bin_global, self)
        return self.rand1bin_global
