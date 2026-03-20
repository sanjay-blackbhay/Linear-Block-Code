# Linear-Block-Code
# Aim
Write a simple python program to Generate Matrix, Codeword, Hamming weight, Syndrome matrix and find the error on received codeword using Linear block code. 
# Tools required
PC & Google Colab
# Theory
A Linear Block Code is an error-correcting code in which k data bits are encoded into n-bit codewords using linear combinations, so that errors during transmission can be detected and corrected.Each block of k bits is converted into an n-bit codeword using a generator matrix. The codewords follow linear properties (sum of any two codewords is also a valid codeword).
# Program
```python
import numpy as np
pb = [] 
Ik = [] 
p = []
m = []
h = []
h_dis = []
r_code = []
err = []
col = int(input("Enter the Parity bits : "))
row = int(input("Enter the Message bits : "))

for i in range (row):
    p = list(map(int, input(f"Enter the row values : {i+1} (Separated by space) : ").split()))  
    pb.append(p)
p_mat = np.array(pb, dtype=int)
Ik=np.eye(row, dtype=int) 
g_mat = np.hstack((p_mat,Ik)) 
n, k = g_mat.T.shape

m = np.array([[1 if (i >> (k - j - 1)) & 1 else 0 for j in range(k)] for i in range(2**k)])

c = np.mod(np.dot(m, g_mat), 2)
for i, row in enumerate(c):
    h_dis1 = np.sum(row)  
    h_dis.append(h_dis1)
h_mat = np.array(h_dis).reshape(1,-1)

d_min = np.min(np.sum(c[1:], axis=1))

h = p_mat[:, :3]
hp = np.hstack((np.eye(n-k, dtype=int), h.T))
ht = hp.T
print('**********')
print('The Generator Matrix is: ')

for r in g_mat: 
    print(" ".join(map(str, r)))
print('**********')
print(f'Message Bits  Codeword   Hamming Weight')
code_word = np.hstack((m, c, h_mat.T))
for r in range(code_word.shape[0]):
    format_row = " ".join(map(str, code_word[r, :k])) + '\t' + " ".join(map(str, code_word[r, k:n+k])) + '\t' + str(code_word[r, -1])
    print(format_row)
print('**********')
print(f'Minimum Hamming distance : {d_min}')

print('**********')
print(f'Parity Check Matrix')
for r in hp:
    print(" ".join(map(str, r)))
print('**********')
print(f'Parity Check Matrix Transpose')
for r in ht:
    print(" ".join(map(str, r)))

rc = list(map(int, input(f"Enter the error codeword : ").split()))  
r_code.append(rc)
r_c = np.array(r_code)

e = np.mod(np.dot(r_c, ht), 2)

print('**********')
print(f"Syndeome of given received codeword is : " + " ".join(map(str, e[0])))
print('**********')
print(f'Syndrome Matrix')
for i in range(n):
    combined_row = np.concatenate((ht[i, :], np.eye(n, dtype=int)[i,:]))
    formatted_row = " ".join(map(str, combined_row[:3])) + '\t' + " ".join(map(str, combined_row[k:]))
    print(f'{formatted_row}')

for i in range(n):
    if np.array_equal(e[0], ht[i, :]):
        err = np.eye(n, dtype=int)[i,:]
print(f"The error postion is : " + " ".join(map(str, err)))

add = err + rc
print(f"The correct codeword is : " + " " .join(map(str,add)))
```
# Output 

<img width="337" height="588" alt="Screenshot 2026-03-13 134156" src="https://github.com/user-attachments/assets/7b2e21d1-39a0-4b8a-8e90-5226bfb98823" />

# Results
Thus the Linear Block Code is verified using Python Program in Google Colab 
