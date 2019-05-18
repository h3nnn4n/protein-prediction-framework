import math
import random
import numpy as np


class LocalitySensitiveHashing:
    def __init__(self, de=None):
        self.de = de

        self.do_lsh = False
        self.n_hashes = None
        self.n_buckets = None
        self.hashes = None
        self.hash_values = None
        self.active_hash1 = 0
        self.active_hash2 = 1
        self.tmp1 = None
        self.tmp2 = None
        self.update_interval = None
        self.change_interval = None

        self.pop_size = self.de.pop_size

    def inject_parameters(self):
        self.do_lsh = self.de.do_lsh
        self.n_hashes = self.de.n_hashes
        self.n_buckets = self.de.n_buckets
        self.update_interval = self.de.update_interval
        self.change_interval = self.de.change_interval

    def lsh_step(self):
        if self.do_lsh and self.de.it % self.change_interval == 0:
            self.change_hash()

    def create_hashs(self):
        self.hashes = [np.random.randint(100, size=self.de.pop[0].total_number_of_angles) for _ in range(self.n_hashes)]

    def apply_hash(self, debug=False):
        # debug = True

        if self.hash_values is None:
            self.hash_values = [[] for _ in range((self.n_buckets + 1) ** 2 + (self.n_buckets + 1))]
        else:
            for k, v in enumerate(self.hash_values):
                if len(v) > 0:
                    self.hash_values[k] = []

        minh1 = None
        maxh1 = None
        minh2 = None
        maxh2 = None

        if self.tmp1 is None:
            self.tmp1 = np.empty(self.pop_size)

        if self.tmp2 is None:
            self.tmp2 = np.empty(self.pop_size)

        tmp1 = self.tmp1
        tmp2 = self.tmp2

        for i in range(self.pop_size):
            h1 = np.dot(self.hashes[self.active_hash1], self.de.pop[i].angles)
            h2 = np.dot(self.hashes[self.active_hash2], self.de.pop[i].angles)
            tmp1[i] = h1
            tmp2[i] = h2

            if minh1 is None or h1 < minh1:
                minh1 = h1

            if maxh1 is None or h1 > maxh1:
                maxh1 = h1

            if minh2 is None or h2 < minh2:
                minh2 = h2

            if maxh2 is None or h2 > maxh2:
                maxh2 = h2

        r1 = (maxh1 - minh1) / self.n_buckets
        b1 = random.random() * r1

        r2 = (maxh2 - minh2) / self.n_buckets
        b2 = random.random() * r2

        if debug:
            print("r1: %8.3f  b1: %8.3f  minh1: %8.3f  maxh1: %8.3f" % (r1, b1, minh1, maxh1))
            print("r2: %8.3f  b2: %8.3f  minh2: %8.3f  maxh2: %8.3f" % (r2, b2, minh2, maxh2))

        for i in range(self.pop_size):
            a, b = math.floor((tmp1[i] - minh1 + b1) / r1), math.floor((tmp2[i] - minh2 + b2) / r2)
            v = (self.n_buckets) * a + b

            try:
                self.hash_values[v].append(i)
            except Exception:
                raise IndexError(v, len(self.hash_values))

        for n, i in enumerate(self.hash_values):
            if debug:
                if len(i) > 0:
                    print(n, i)

        if debug:
            print()

    def change_hash(self):
        old1 = self.active_hash1
        old2 = self.active_hash2

        while old1 == self.active_hash1:
            self.active_hash1 = random.randint(0, self.n_hashes - 1)

        while old2 == self.active_hash2 and self.active_hash1 == self.active_hash2:
            self.active_hash2 = random.randint(0, self.n_hashes - 1)
