#include "output.hpp"

void FileOutput::writeDemesHeader() {
    file << "Generation,Deme,Side,Population,OriginTime,AverageArray" << std::endl;
}
void FileOutput::writeDemesFile(Tumour& tumour) {
    for (int i = 0; i < tumour.getNumDemes(); i++) {
      file << tumour.getGensElapsed() << "," << i << ","
           << tumour.getDeme(i).getSide() << ","
           << tumour.getDeme(i).getPopulation() << ","
           << tumour.getDeme(i).getOriginTime() << ",";
      for (int j = 0; j < tumour.getDeme(i).getAverageArray().size(); j++) {
        file << tumour.getDeme(i).getAverageArray()[j] << ";";
        }
        file << std::endl;
    }
}

void FileOutput::writeCellsFile(Tumour& tumour) {
    return;
}
