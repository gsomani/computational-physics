import random

def dice_roll():	
	return random.randint(1,6),random.randint(1,6) 

n = 2**20
count = 0
uniform_prob = 1/36

for i in range(n):
	d1,d2 = dice_roll()
	count += d1==6 and d2==6

prob = count/n
diff = abs(count - uniform_prob) 
print("Probabilty (by Monte carlo) = %.3f" % prob)
print("Uniform distribution probability = %.3f" % uniform_prob)
print("Difference = %.3f" % abs(uniform_prob-prob))
