import numpy as np
import pandas as pd
import math
import sys
import mysymnmf

''' from now on-
W = normalized similarity matrix (laplacian) (1.3)
A = The Similarity Matrix (1.1)
D = The diagonal degree Matrix (1.2)
M = random matrix sign
H = decomposition matrix'''

np.random.seed(0)


# norm calculators----------------------------------------
def F_norm(M):
    # M is numpy array
    ans = np.sum(M**2)
    return math.sqrt(ans)

def euclidian_dist(x1, x2):
    ans = 0
    for i in range(len(x1)):
        ans += (x1[i]-x2[i])**2
    return ans


# matrix builders-------------------------------------------

def build_similarity_mat_A(points):
    # build A from points array
    n = len(points)
    A = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i != j:
                A[i,j] = math.exp(-0.5 * euclidian_dist(points[i,:], points[j,:]))
    return A


def build_diag_deg_D(A):
    # build D from A
    v = np.sum(A, axis=1)
    return np.diagflat([v])


def build_laplacian_W(A, D):
    # A = Similarity matrix, D = diagonal degree mat
    # build W
    D_new = np.zeros(D.shape)
    for i in range(len(D)):
        D_new[i,i] = D[i,i] ** -0.5
    tmp = np.matmul(D_new, A)
    W = np.matmul(tmp, D_new)
    return W

def init_H(W, k):
    # initialize H
    # W = laplasian matrix, k = num of clusters
    m = np.mean(W)
    H = np.random.uniform(low=0, high=2*math.sqrt(m/k), size=(len(W),k))
    return H

def one_update_H(H_i, W, beta=0.5):
    # perform one iteration of updating H
    a = np.matmul(W,H_i)
    b = np.matmul(H_i, H_i.T)
    b = np.matmul(b,H_i)
    new = a/b
    new *= beta
    new += (1-beta)
    return H_i * new

def calculate_H(H, W, iter, epsilon = 0.01):
    # update H until convergance / max iter is reached
    H_new = one_update_H(H, W)
    norm = F_norm((H_new - H))
    while norm >= epsilon or iter > 0:
        H_old = H_new
        H_new = one_update_H(H_old, W)
        norm = F_norm((H_new - H_old))
        iter -= 1
    return H_new


# main functions-------------------------------------------------------------
def run_python(goal, points, k, it=300, eps=0.01):
    # Perform the symNMF algorithm according to the goal - with python functions - test
    A = build_similarity_mat_A(points)
    if goal == "sym":
        return A
    D = build_diag_deg_D(A)
    if goal == "ddg":
        return D
    W = build_laplacian_W(A,D)
    if goal == "norm":
        return W
    H_0 = init_H(W, k)
    H = calculate_H(H_0, W, it, eps)
    return H

def run(goal, points, k):
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
    

def main():
    try:
        k = int(sys.argv[1])
        goal = sys.argv[2]
        file_path = sys.argv[3]
        df = pd.read_csv(file_path, sep=",", header=None)
        points = df.values
        n = points.shape[0]
        d = points.shape[1]
        H = run(goal, points.tolist(), k)
        for i in range(len(H)):
            s = [f'{num:.4f}' for num in H[i]]
            print(*s, sep = ",") 

    except Exception as e:
        print("An Error Has Occurred")
        print(e)
        exit(1)



if __name__ == '__main__':
    main()
