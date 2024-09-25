CC = gcc

CFLAGS = -Wall -g

TARGET = pnorm

all: $(TARGET)

$(TARGET): pnorm.o
	$(CC) $(CFLAGS) -o $(TARGET) pnorm.o

pnorm.o: pnorm.c
	$(CC) $(CFLAGS) -c pnorm.c

clean:
	rm -f *.o $(TARGET)

edit:
	vim pnorm.c
