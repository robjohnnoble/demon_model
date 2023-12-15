#include "output.hpp"

void FileOutput::write(Tumour& tumour) {
    file << "Deme\tSide\tPopulation\tAverageArray" << std::endl;
    for (int i = 0; i < tumour.getNumDemes(); i++) {
        file << i << "\t" << tumour.getDeme(i).getSide() << tumour.getDeme(i).getPopulation() << "\t";
        for (int j = 0; j < tumour.getDeme(i).getAverageArray().size(); j++) {
            file << tumour.getDeme(i).getAverageArray()[j] << " ";
        }
        file << std::endl;
    }
}