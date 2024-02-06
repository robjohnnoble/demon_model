# Compiler settings for debugging
CXX = g++
CXXFLAGS = -Wall -g -O0 -std=c++11 -I include/methdemon -I /opt/homebrew/Cellar/boost/1.84.0/include/
LDFLAGS =

# Directories
SRCDIR = src
INCDIR = include/methdemon
BINDIR = bin
LOGDIR = logs

# Find all source files in the source directory
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
# Replace .cpp from SOURCES with .o to get object files
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(BINDIR)/%.o)

# Name of the executable
EXECUTABLE = $(BINDIR)/methdemon

# Makefile targets
all: $(LOGDIR) $(BINDIR) $(EXECUTABLE)

$(LOGDIR):
	mkdir -p $(LOGDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^ 2>&1 | tee -a $(LOGDIR)/linking.log

$(BINDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ 2>&1 | tee -a $(LOGDIR)/$(notdir $<).log

clean:
	rm -rf $(BINDIR) $(LOGDIR)

