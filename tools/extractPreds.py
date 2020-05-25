preds = []

with open("twz_5D.out", "r") as a_file:
    for line in a_file:
        if("Cross-section" in line):
            token = line.split(" ")
            pred = []
            pred.append(float(token[9]))
            preds.append(pred)
        
print("predictions = " )
print(preds)
