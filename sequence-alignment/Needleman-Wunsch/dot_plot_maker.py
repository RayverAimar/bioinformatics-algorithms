import numpy as np
from PIL import Image

NUMBER_OF_SEQUENCES = 3

def read_txt():
    data = []
    path = "processed_sequences.txt"
    try:
        with open(path, 'r') as file:
            for i in range(NUMBER_OF_SEQUENCES):
                name = next(file)
                sequence = next(file)
                data.append((name.strip(), sequence.strip()))
        return data
    except FileNotFoundError:
        print(f"The file '{path}' was not found.")
        return data

def get_dot_plots(sequences):
    def get_dot_plot(tuple_a, tuple_b):
        sequence_a = tuple_a[1]
        sequence_b = tuple_b[1]
        image_array = np.zeros((len(sequence_a), len(sequence_b)))
        for i in range(image_array.shape[0]):
            for j in range(image_array.shape[1]):
                if sequence_a[i] != sequence_b[j]:
                    image_array[i, j] = 255
        image = Image.fromarray(image_array)
        image = image.convert("RGB")
        image_path = "./dot_plots/" + tuple_a[0] + '_x_' + tuple_b[0] + '.jpg'
        image.save(image_path)
        print("Succesfully saved image in:", image_path)
    for i in range(len(sequences)):
        for j in range(i, len(sequences)):
            get_dot_plot(sequences[i], sequences[j])

sequences = read_txt()
get_dot_plots(sequences)