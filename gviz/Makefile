CC = g++
CFLAGS = -m64 -std=c++11
LIB = -lGL -lGLU -lglfw3 -lX11 -lXxf86vm -lXrandr -lpthread -lXi

NAME = gviz

default: all
all: clean $(NAME)

$(NAME).o: $(NAME).cpp
	$(CC) $(CFLAGS) -c $(NAME).cpp
	
$(NAME): $(NAME).o
	$(CC) $(NAME).o -o $(NAME) $(LIB)

clean:
	rm -f $(NAME) $(NAME).o