import numpy as np

A = "ACTGTGCGTACGTGA"
B = "AGCTAGTATACGGAT"

# A = "TGGTG"
# B = "ATCGT"

# A = "CGTGAATTCAT"
# B = "GACTTAC"

F = np.zeros((len(A)+1,len(B)+1), dtype=int)
print(F)
# print("np.zeros((5,5)) :      ", F)
# print("Nombre de dimensions : ", F.ndim)  # Nombre de dimensions 2
# print("Nombre d'éléments :    ", F.size)  # Nombre d'éléments 5*5
# print("Shape :                ", F.shape) # Shape : (5,5)
# print("np.arange(6) :         ", np.arange(6))

# Gap penalty

# d = 

# for i in range(0,len(A)):
#     pass

match = 1
mismatch = -1
gap = -1

for i in range(0, len(A)+1):
    F[i,0] = gap * i
for j in range(0, len(B)+1):
    F[0,j] = gap * j

for i in range(1, len(A)+1):
    for j in range(1, len(B)+1):
        test = F[i-1,j-1] + (match if A[i-1] == B[j-1] else mismatch)
        left = F[i-1,j] + gap
        up = F[i,j-1] + gap
        F[i,j] = max(test, left, up)

# F = np.flip(F)
print(np.flip(F))



AlignmentA = ""
AlignmentB = ""
i = len(A)-1
j = len(B)-1
while i >= 0 and j >= 0:
    if i > 0 and j > 0 and F[i, j] == F[i-1, j-1] + (match if A[i] == B[j] else mismatch):
        AlignmentA = A[i] + AlignmentA
        AlignmentB = B[j] + AlignmentB
        i = i - 1
        j = j - 1
    elif i > 0 and F[i, j] == F[i-1, j] + gap:
        AlignmentA = A[i] + AlignmentA
        AlignmentB = "-" + AlignmentB
        i = i - 1
    else:
        AlignmentA = "-" + AlignmentA
        AlignmentB = B[j] + AlignmentB
        j = j - 1

print(AlignmentA)
print(AlignmentB)


# def NWScore(X, Y):
#     Score(0, 0) = 0 # 2 * (length(Y) + 1) array
#     for j = 1 to length(Y)
#         Score(0, j) = Score(0, j - 1) + Ins(Yj)
#     for i = 1 to length(X) # Init array
#         Score(1, 0) = Score(0, 0) + Del(Xi)
#         for j = 1 to length(Y)
#             scoreSub = Score(0, j - 1) + Sub(Xi, Yj)
#             scoreDel = Score[0, j] + Del(Xi)
#             scoreIns = Score[1, j - 1] + Ins(Yj)
#             Score[1, j] = max(scoreSub, scoreDel, scoreIns)
#         # Copy Score[1] to Score[0]
#     for iitem, item in enumerate(Score[1]):
#         Score[0, iitem] = item
#     for j in range(len(Y)):
#         LastLine[j] = Score[1, j]
#     return LastLine



# def Hirschberg(X, Y):
#     Z = ""
#     W = ""
#     if len(X) == 0:
#         for i in range(1, len(Y)):
#             Z = Z + '-'
#             W = W + Y[i]
#     elif len(Y) == 0:
#         for i in range(1, len(X)):
#             Z = Z + X[i]
#             W = W + '-'
#     elif len(X) == 1 or len(Y) == 1:
#         (Z, W) = NeedlemanWunsch(X, Y)
#     else
#         xlen = length(X)
#         xmid = length(X) / 2
#         ylen = length(Y)

#         ScoreL = NWScore(X1:xmid, Y)
#         ScoreR = NWScore(rev(Xxmid+1:xlen), rev(Y))
#         ymid = arg max ScoreL + rev(ScoreR)

#         (Z,W) = Hirschberg(X1:xmid, y1:ymid) + Hirschberg(Xxmid+1:xlen, Yymid+1:ylen)
#     end
#     return (Z, W)

