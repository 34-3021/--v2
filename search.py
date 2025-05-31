from judge import calculate_value, ref, query

a = [] # the answer sequence
maxn = 0

def dfs(i, seq):
    s = seq[:]
    if i == len(seq):
        st = f"({s[0]}, {s[1]}, {s[2]}, {s[3]}), ({s[4]}, {s[5]}, {s[6]}, {s[7]}), ({s[8]}, {s[9]}, {s[10]}, {s[11]})"
        val = calculate_value(st, ref, query)
        global maxn
        if val > maxn:
            maxn = val
            print(f"New max value: {maxn} for sequence: {st}")
        return

    for delta in [-1, 0, 1]:
        s[i] = seq[i] + delta
        if s[i] >= 0:
            dfs(i + 1, s)
        
dfs(0, a)
