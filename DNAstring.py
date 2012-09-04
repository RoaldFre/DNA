import random
alphabet = 'ACTG'
total = 100
string=''

for x in range(100):
    i=random.sample(alphabet,1)
    i=str(i).strip("'[]'")
    string+=i
    
print(string)
