#include "output.hpp"

void FileOutput::writeDemesHeader() {
    file << "Generation\tDeme\tSide\tPopulation\tOriginTime\tAverageArray" << std::endl;
}
void FileOutput::writeDemesFile(Tumour& tumour) {
    for (int i = 0; i < tumour.getNumDemes(); i++) {
        file << tumour.getGensElapsed() << "\t" << i << "\t" << tumour.getDeme(i).getSide() << "\t" << tumour.getDeme(i).getPopulation() << "\t" << tumour.getDeme(i).getOriginTime() << "\t";
        for (int j = 0; j < tumour.getDeme(i).getAverageArray().size(); j++) {
            file << tumour.getDeme(i).getAverageArray()[j] << " ";
        }
        file << std::endl;
    }
}

void FileOutput::writeCellsFile(Tumour& tumour) {
    return;
}
