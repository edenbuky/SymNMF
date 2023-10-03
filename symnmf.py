import numpy as np
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



def init_H(W, k):
    # initialize H
    # W = laplasian matrix, k = num of clusters
    m = np.mean(W)
    H = np.random.uniform(low=0, high=2*math.sqrt(m/k), size=(W.shape[0],k))
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
        points = []
        with open(file_path, 'r') as file:
            for line in file:
                point = [float(num) for num in line.split(',')]
                points.append(point)
        H = run(goal, points, k)
        for i in range(len(H)):
            s = [f'{num:.4f}' for num in H[i]]
            print(*s, sep = ",") 

    except Exception as e:
        print("An Error Has Occurred")
        print(e)
        exit(1)



if __name__ == '__main__':
    main()
