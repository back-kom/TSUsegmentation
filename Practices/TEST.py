# S = [1,2,6,5,5,5,3]
# E = [5,5,7,6,6,7,8]
# n = len(S)
# S.sort()
# E.sort()
# count = 1
# j=0
# for i in range(1, n):
#     if j < i and S[i] >= E[j]:
#         j += 1
#     elif S[i] < E[j]:
#         count += 1
# print(count)

# import math
#
#
# class Solution:
#
#     k=4
#     x=[-1,2,-4,2,4]
#     y=[1,2,-4,2,1]
#     comb = list(zip(x,y))
#     s=[]
#     for i,j in zip(x,y):
#         r = pow(i,2)+pow(j,2)
#         s.append(r)
#     s.sort()
#     print(math.ceil(math.sqrt(s[k-1])))


A=[5,1,3,4,2]
n=len(A)
next_odd=[0]*n
next_even=[0]*n
stack = []
index_sorted = sorted(range(n), key = lambda i:A[i])









