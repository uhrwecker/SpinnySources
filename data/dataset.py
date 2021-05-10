from data.trajectory import TrajectoryResults
import os
import time

class Dataset:
    def __init__(self, s, path, both=True):
        self.s = s
        self.path = path

        self.all = None

        start = time.time()
        if both:
            print(f'Start reading data for s = {self.s} ...')
            print(f'- Read from 0          to pi / 2     ...')
            self.one_up, self.one_down = self.read_datapoints(path+'1/')

            print(f'- Read from pi / 2     to pi         ...')
            self.two_up, self.two_down = self.read_datapoints(path + '2/')

            print(f'- Read from pi         to 3 * pi / 2 ...')
            self.three_up, self.three_down = self.read_datapoints(path + '3/')

            print(f'- Read from 3 * pi / 2 to 2 * pi     ...')
            self.four_up, self.four_down = self.read_datapoints(path + '4/')
        else:
            self.all = self.read_datapoints(path)

        duration = time.time() - start
        print(f'Done! Took {duration}s (or {duration / 60}m).')

    def read_datapoints(self, path):
        up = {}
        down = {}
        for file in os.listdir(path):
            if file.endswith('.ini'):
                continue
            fp = file[:-4]

            point = TrajectoryResults(path + fp)

            if fp.startswith('up'):
                up[str(point.phi0)] = point
            elif fp.startswith('down'):
                down[str(point.phi0)] = point
            else:
                up[fp] = point

        return up, down

    def get_all_data(self):
        return self.one_up, self.one_down, self.two_up, self.two_down, self.three_up, self.three_down, self.four_up, self.four_down

    def get_all_up(self):
        return self.one_up, self.two_up, self.three_up, self.four_up

    def get_all_down(self):
        return self.one_down, self.two_down, self.three_down, self.four_down

    def get_all(self):
        return self.all
