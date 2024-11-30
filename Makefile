# ~~ Compiler and target ~~
CXX := g++
TARGET := main


# ~~ Directories ~~
SRC_DIR := src
INC_DIR := inc
OBJ_DIR := bin

# ~~ ROOT flags and compiler flags ~~
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
CXXFLAGS := -Wall -I$(INC_DIR) $(ROOTCFLAGS)

# ~~ Source and object files ~~
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRCS))
MAIN_SRC := $(TARGET).cpp
MAIN_OBJ := $(OBJ_DIR)/$(TARGET).o

# === Build Target ===
$(TARGET): $(MAIN_OBJ) $(OBJS)
	$(CXX) -o $@ $^ $(ROOTLIBS)

# === Compile Main Program ===
$(MAIN_OBJ): $(MAIN_SRC) | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# === Compile Modules ===
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ~~ Ensure bin directory exists ~~
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# === Clean Up ===
.PHONY: clean
clean:
	rm -rf $(OBJ_DIR) $(TARGET)
