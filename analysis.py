from sklearn.metrics import silhouette_score
import numpy as np
import math
import sys
import pandas as pd
import mysymnmf

epsilon = 0.0001
iter = 300

def init_H(W, k):
    # initialize H
    # W = laplasian matrix, k = num of clusters
    m = np.mean(W)
    H = np.random.uniform(low=0, high=2*math.sqrt(m/k), size=(len(W),k))
    return H


def get_clusters(H):
    #H list of lists
    H = np.array(H)
    clusters = np.argmax(H, axis=1) # clusters[i] = num of cluster point i is in
    




def distance(p1, p2):
    sum = 0
    for i in range(len(p1)):
        tmp = p1[i] - p2[i]
        sum += (tmp**2)
    return sum ** 0.5

class K_Means:
    def __init__(self, k, iter, points):
        self.k = k
        self.iter = iter
        self.points = points
        self.centroids = [c for c in points[:k]]
        self.clusters = [[c] for c in self.centroids]
        self.mk = [False for i in range(k)]
        self.d = len(self.centroids[0])


    def calc_mean(self, cluster):
        new = [0 for x in range(self.d)]
        for point in cluster:
            for i in range(len(new)):
                new[i] += point[i]
        return [new[j]/len(cluster) for j in range(len(new))]


    def run_iter(self):
        global epsilon
        # assign points to clusters
        for i in range(len(self.points)):
            point = self.points[i]
            dis = float('inf')
            idx = 0
            for cluster_idx in range(len(self.centroids)):
                cluster = self.centroids[cluster_idx]
                tmp = distance(point, cluster)
                if tmp < dis:
                    dis = tmp
                    idx = cluster_idx
            self.clusters[idx].append(point)
        # update centroids
        for j in range(len(self.centroids)):
            new_mean = self.calc_mean(self.clusters[j])
            dmk = distance(new_mean, self.centroids[j])
            if dmk < epsilon:
                self.mk[j] = True
            self.centroids[j] = new_mean
        # erase points from clusters
        self.clusters = [[] for i in range(self.k)]

    def run(self, pnt):
        iters = self.iter
        while iters > 0 and not all(self.mk):
            self.run_iter()
            iters -= 1
        if pnt:
            for centroid in self.centroids:
                print(','.join(['%.4f' % coord for coord in centroid]))



def main():
    try:
        k = int(sys.argv[1])
        file_path = sys.argv[2]
        df = pd.read_csv(file_path, sep=",", header=None)
        points = df.values
        n = points.shape[0]
        W = mysymnmf.norm(points)
        H0 = init_H(np.array(W), k)
        H = mysymnmf.symnmf(points, H0.tolist(), W, n, k)

        kmeans = K_Means(k, iter, points)
        kmeans.run()

        kmeans_clusters = kmeans.centroids

    except Exception as e:
        print("An Error Has Occurred")
        # print(e)
        exit(1)



if __name__ == '__main__':
    main()
