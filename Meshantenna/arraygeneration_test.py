import numpy as np
import random
from collections import deque

def generate_array_with_min_path(params, ARRAY_SHAPE=(10,10), min_path_length=5):
    """
    Generate array_2d so that the path from the 'ant' in row -1
    to any 1 in row -2 is at least `min_path_length` steps if it exists.
    If a path < min_path_length is found, we regenerate.
    """
    while True:
        array_2d = np.ones(ARRAY_SHAPE, dtype=int)

        # 1. Randomly fill 0/1 with p=0.4 for "tiles"
        for i in range(ARRAY_SHAPE[0]):
            for j in range(ARRAY_SHAPE[1]):
                array_2d[i][j] = 1 if (random.random() < 0.4) else 0

        # 2. Randomly fill second-to-last row with p=0.1
        #    (This is row index -2 in Python, but let's be explicit: row = ARRAY_SHAPE[0] - 2)
        for col in range(ARRAY_SHAPE[1]):
            array_2d[ARRAY_SHAPE[0]-2][col] = 1 if (random.random() < 0.1) else 0

        # 3. Set the last row to 0, then place the "ant" = 1 in the specified column
        array_2d[-1] = 0
        if params['ant_fp'] == -1:
            # If random
            ant_col = random.randint(0, ARRAY_SHAPE[1] - 1)
        else:
            ant_col = params['ant_fp']
        array_2d[-1][ant_col] = 1

        # 4. Check BFS distance from (last_row, ant_col) to any "1" in row -2
        start = (ARRAY_SHAPE[0]-1, ant_col)
        distance = find_min_distance(array_2d, start, target_row=ARRAY_SHAPE[0]-2)
        
        # If there is NO path, distance will be None (or we can interpret as infinite)
        # If there is a path, we want it to be >= min_path_length
        if distance is None or distance >= min_path_length:
            # Good scenario: either no path or path >= 5
            return array_2d
        else:
            # We found a path < min_path_length, so re-generate
            continue


def find_min_distance(array_2d, start, target_row):
    """
    BFS to find the shortest distance from `start`
    to any 1 in `target_row`.
    
    array_2d: 2D numpy array of 0/1
    start: (row, col) coordinate to start BFS
    target_row: integer row index for the row
    that we want to reach (where array_2d[row][col] == 1).
    
    Return the minimum distance if found, else None.
    Distance is counted in "steps" on a 4-connected grid.
    """
    nrows, ncols = array_2d.shape
    sr, sc = start
    
    # If start itself is not walkable, or out of range, no path
    if not (0 <= sr < nrows and 0 <= sc < ncols and array_2d[sr][sc] == 1):
        return None

    visited = np.zeros_like(array_2d, dtype=bool)
    visited[sr][sc] = True
    
    # BFS queue holds tuples: (row, col, distance)
    queue = deque([(sr, sc, 0)])
    
    # Directions for up, down, left, right
    directions = [(0,1), (0,-1), (1,0), (-1,0)]
    
    while queue:
        r, c, dist = queue.popleft()
        
        # Check if we've reached the target row at a cell that is 1
        if r == target_row and array_2d[r][c] == 1:
            return dist
        
        # Otherwise, expand neighbors
        for dr, dc in directions:
            rr, cc = r + dr, c + dc
            if 0 <= rr < nrows and 0 <= cc < ncols:
                if not visited[rr][cc] and array_2d[rr][cc] == 1:
                    visited[rr][cc] = True
                    queue.append((rr, cc, dist + 1))
    
    # If we exhaust the queue without reaching target row, no path
    return None


if __name__ == "__main__":
    params = {
        'ant_fp': -1  # or set some integer column
    }
    from time import time
    now = time()
    ARRAY_SHAPE = (10, 10)
    array_2d = generate_array_with_min_path(params, ARRAY_SHAPE, min_path_length=5)
    print(array_2d)
    print("Time taken:", time() - now)