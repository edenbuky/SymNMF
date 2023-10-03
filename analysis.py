from sklearn.metrics import silhouette_score
import numpy as np
import math
import sys
import mysymnmf


def init_H(W, k):
    # initialize H
    # W = laplasian matrix, k = num of clusters
    m = np.mean(W)
    H = np.random.uniform(low=0, high=2*math.sqrt(m/k), size=(len(W),k))
    return H

def run_symnmf(goal, points, k):
    # Perform the symNMF algorithm according to the goal
    if goal == "sym":
        return mysymnmf.sym(points)
    elif goal == "ddg":
        return mysymnmf.ddg(points)
    elif goal == "norm":
        return mysymnmf.norm(points)
    else:
        W = mysymnmf.norm(points)
        H = init_H(np.array(W),k)
        n = H.shape[0]
        H_new = H.tolist()
        return mysymnmf.symnmf(n, k, H_new, W)


def get_clusters_symnmf(goal, points, k):
    # in: points as list of lists, out: np array
    H = run_symnmf(goal, points, k)
    H = np.array(H)
    clusters_idx = np.argmax(H, axis=1) # clusters[i] = num of cluster point i is in
    return clusters_idx

    
def get_clusters_kmeans(k, points):
    # in: points as list of lists, out: np array
    kmeans = K_Means(k, points)
    kmeans.run_kmeans()
    return np.array(kmeans.cluster_labels)


def distance(p1, p2):
    sum = 0
    for i in range(len(p1)):
        tmp = p1[i] - p2[i]
        sum += (tmp**2)
    return sum ** 0.5

class K_Means:
    def __init__(self, k, points, iter=300, epsilon=0.0001):
        self.k = k
        self.iter = iter
        self.epsilon = epsilon
        self.points = points
        self.cluster_labels = [i for i in range(len(points))]
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
            self.cluster_labels[i] = idx
        # update centroids
        for j in range(len(self.centroids)):
            new_mean = self.calc_mean(self.clusters[j])
            dmk = distance(new_mean, self.centroids[j])
            if dmk < self.epsilon:
                self.mk[j] = True
            self.centroids[j] = new_mean
        # erase points from clusters
        self.clusters = [[] for i in range(self.k)]

    def run_kmeans(self, pnt=False):
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
        goal = "symnmf"
        points = []
        with open(file_path, 'r') as file:
            for line in file:
                point = [float(num) for num in line.split(',')]
                points.append(point)

        symnmf_labels = get_clusters_symnmf(goal, points, k)
        kmeans_labels = get_clusters_kmeans(k, points)

        points = np.array(points)

        symnmf_score = silhouette_score(points, symnmf_labels)
        kmeans_score = silhouette_score(points, kmeans_labels)

        print("nmf: %.4f" % symnmf_score)
        print("kmeans: %.4f" % kmeans_score)



    except Exception as e:
        print("An Error Has Occurred")
        # print(e)
        exit(1)



if __name__ == '__main__':
    main()
