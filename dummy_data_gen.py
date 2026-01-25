import csv
import random

def generate_recidivism_data(filename, num_rows):
    headers = ['week', 'arrest', 'fin', 'age', 'race', 'wexp', 'mar', 'paro', 'prio']
    
    print(f"generating...")

    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(headers)

        for _ in range(num_rows):
            week = random.randint(1, 52)          # Week of arrest/observation
            arrest = random.choice([0, 1])        # 1 if arrested, 0 otherwise
            fin = random.choice([0, 1])           # Financial aid received
            age = random.randint(18, 45)          # Realistic age range
            race = random.choice([0, 1])          # Binary race indicator
            wexp = random.choice([0, 1])          # Prior work experience
            mar = random.choice([0, 1])           # Married (1) or not (0)
            paro = random.choice([0, 1])          # On parole
            prio = random.randint(0, 15)          # Number of prior convictions

            writer.writerow([week, arrest, fin, age, race, wexp, mar, paro, prio])

    print(f"Success! Data saved to {filename}")

output_file = 'dummy_data.csv'
total_rows = 10_000_000

generate_recidivism_data(output_file, total_rows)
