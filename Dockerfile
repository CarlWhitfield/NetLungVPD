FROM ubuntu:latest

# Install necessary dependencies, including MPI libraries
RUN apt-get update && apt-get install -y \
    git \
    libeigen3-dev \
    g++

# Set the working directory inside the container
WORKDIR /app
# Copy the source code and necessary files to the container
RUN mkdir /app/Tube-network-headers
RUN git clone https://github.com/CarlWhitfield/Tube-network-headers /app/Tube-network-headers
COPY ./include /app/include
COPY ./src /app/src

RUN ls -a /app

# Build the C++ program using g++
RUN g++ -o OOLungSim ./src/*.cpp -I./include/ -I./Tube-network-headers -I/usr/include/eigen3 -I/usr/include/boost -std=c++11 -O3

# Set the entry point for the container
ENTRYPOINT ["/app/OOLungSim"]
