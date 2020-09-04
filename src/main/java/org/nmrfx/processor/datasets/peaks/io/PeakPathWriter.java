/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.nmrfx.processor.datasets.peaks.io;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import org.nmrfx.processor.datasets.peaks.InvalidPeakException;
import org.nmrfx.processor.datasets.peaks.PeakList;
import org.nmrfx.processor.datasets.peaks.PeakPath;
import org.nmrfx.processor.datasets.peaks.PeakPath.PeakDistance;

/**
 *
 * @author brucejohnson
 */
public class PeakPathWriter {

    public void writeToSTAR3(FileWriter chan, PeakPath peakPath, int id) throws IOException, InvalidPeakException {
        char stringQuote = '"';
        List<PeakList> peakLists = peakPath.getPeakLists();
        PeakList firstList = peakLists.get(0);

        chan.write("save_" + peakPath.getName() + "\n");
        chan.write("_NMRFx_peak_path.Sf_category                 ");
        chan.write("nmrfx_peak_path\n");
        chan.write("_NMRFx_peak_path.Sf_framecode                 ");
        chan.write(peakPath.getName() + "\n");
        chan.write("_NMRFx_peak_path.ID                          ");
        chan.write(id + "\n");
        chan.write("_NMRFx_peak_path.Type               ");
        chan.write(peakPath.getPathMode().toString().toLowerCase());
        chan.write("\n");
        chan.write("_NMRFx_peak_path.Sample_ID                   ");
        chan.write(".\n");
        int nDim = firstList.nDim;
        chan.write("_NMRFx_peak_path.Number_of_spectral_dimensions ");
        chan.write(String.valueOf(nDim) + "\n");
        chan.write("_NMRFx_peak_path.Units                       ");
        chan.write(peakPath.getUnits() + "\n");
        chan.write("_NMRFx_peak_path.Details                       ");
        if (peakPath.getDetails().length() != 0) {
            chan.write(stringQuote + peakPath.getDetails() + stringQuote + "\n");
        } else {
            chan.write(".\n");
        }
        chan.write("\n");

        List<String> loopStrings = peakPath.getSTAR3LoopStrings();
        chan.write("loop_\n");
        for (String loopString : loopStrings) {
            chan.write(loopString + "\n");
        }
        chan.write("\n");
        int nPeakLists = peakLists.size();
        for (int j = 0; j < nPeakLists; j++) {
            chan.write(peakPath.getSTAR3String(j) + "\n");
        }
        chan.write("stop_\n");
        chan.write("\n");

        String[] dimLoops = {"_Fit_criteria.Spectral_dim_ID", "_Fit_criteria.Spectral_dim_Scale", "_Fit_criteria.Spectral_dim_Tolerance"};
        chan.write("loop_\n");
        for (String loopString : dimLoops) {
            chan.write(loopString + "\n");
        }
        chan.write("\n");

        for (int iDim = 0; iDim < nDim; iDim++) {
            chan.write(peakPath.getSTAR3DimString(iDim));
            chan.write("\n");
        }
        chan.write("stop_\n");
        chan.write("\n");

        loopStrings = peakPath.getSTAR3PathLoopStrings();
        chan.write("loop_\n");
        for (String loopString : loopStrings) {
            chan.write(loopString + "\n");
        }
        chan.write("\n");
        Collection<PeakPath.Path> paths = peakPath.getPaths();
        int elem = 1;

        for (PeakPath.Path path : paths) {
            List<PeakPath.PeakDistance> peakDists = path.getPeakDistances();
            int iList = 0;
            for (PeakDistance peakDist : peakDists) {
                PeakList peakList = peakLists.get(iList++);
                chan.write(String.format("%4d", elem++));
                chan.write(" ");
                chan.write(String.format("%4d", path.getFirstPeak().getIdNum()));
                chan.write(" ");
                chan.write(String.format("%4d", peakList.getId()));
                chan.write(" ");
                if (peakDist != null) {
                    chan.write(String.format("%5d", peakDist.getPeak().getIdNum()));
                } else {
                    chan.write("    ?");
                }
                chan.write("\n");
            }
        }
        chan.write("stop_\n");
        chan.write("\n");

        loopStrings = peakPath.getSTAR3ParLoopStrings();
        chan.write("loop_\n");
        for (String loopString : loopStrings) {
            chan.write(loopString + "\n");
        }
        chan.write("\n");
        elem = 1;
        int useDims = peakPath.getPathMode() == PeakPath.PATHMODE.PRESSURE ? nDim : 1;
        for (PeakPath.Path path : paths) {
            for (int iDim = 0; iDim < useDims; iDim++) {
                chan.write(path.toSTAR3ParString(elem, path.getFirstPeak().getIdNum(), iDim));
                chan.write("\n");
                elem++;
            }
        }
        chan.write("stop_\n");
        chan.write("\n");

        chan.write("\nsave_\n\n");
    }

}
