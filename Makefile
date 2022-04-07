### Compiler & flags
CC=g++
CFLAGS=-c

### Directories
SRC_DIR=src
OBJ_DIR=build
HEAD_DIR=include

### Executable name
EXEC_R2= R2
EXEC_DEF= DEF

### Sources (in src directory)
SOURCES_R2=$(addprefix $(SRC_DIR)/, \
	stt_potential_solver.cpp LSODA.cpp asa047.cpp)

SOURCES_DEF = $(addprefix $(SRC_DIR)/, \
	stt_solver.cpp LSODA.cpp asa047.cpp)

### Objects (in obj directory)
OBJECTS_R2=$(SOURCES_R2:$(SRC_DIR)%.cpp=$(OBJ_DIR)%.o)
OBJECTS_DEF=$(SOURCES_DEF:$(SRC_DIR)%.cpp=$(OBJ_DIR)%.o)

### Rules: #######################################

### General target (executable):
all: $(EXEC_R2) $(EXEC_DEF)

### How to make the executable:
$(EXEC_R2): $(OBJECTS_R2) 
	$(CC) $^ -o $@ -lm

$(EXEC_DEF): $(OBJECTS_DEF) 
	$(CC) $^ -o $@ -lm


### How to make every object:
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) $< -o $@

### How to clean up:
clean: 
	rm $(OBJECTS) $(EXEC_R2) $(EXEC_DEF)
