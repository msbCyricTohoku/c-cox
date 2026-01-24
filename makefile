CC = gcc
CFLAGS = -O3 -march=native -flto -ffast-math -DHAVE_INLINE
LDFLAGS = -lgsl -lgslcblas -lm

TARGET = cox_model
SRC = mpc.c

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET)
