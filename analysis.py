from sklearn.metrics import silhouette_score
import numpy as np
import math

def get_clusters(H):
    # is np array
    clusters = np.argmax(H, axis=1) # clusters[i] = num of cluster point i is in