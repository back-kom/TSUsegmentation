A = [10, 13, 12, 14, 15]
N = len(A)


def make(B):
    ans = [None] * N
    stack = []  # invariant: stack is decreasing
    for i in B:
        while stack and i > stack[-1]:
            ans[stack.pop()] = i
        stack.append(i)
    return ans


B = sorted(range(N), key=lambda i: A[i])
oddnext = make(B)
#B.sort(key=lambda i: -A[i])
B = sorted(range(N), key=lambda i:A[i], reverse=True)
evennext = make(B)

odd = [False] * N
even = [False] * N
odd[N - 1] = even[N - 1] = True

for i in range(N - 2, -1, -1):
    if oddnext[i] is not None:
        odd[i] = even[oddnext[i]]
    if evennext[i] is not None:
        even[i] = odd[evennext[i]]

print(sum(odd))


