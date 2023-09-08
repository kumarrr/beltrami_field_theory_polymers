CC:=		gcc
FLAGS:=		-lm -O2
FILES:=		src/* src/util/*
BIN:=		sphereAND
MAIN:=		src/gauge_sphere_alpha_v5_anderson.c

$(BIN): $(FILES)
	$(CC) $(MAIN) $(FLAGS) -o $@ 
	
clean: 
	rm -rf process_data dat *.dat src/*.c~  *.log
