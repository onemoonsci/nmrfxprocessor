/*
 * NMRFx Processor : A Program for Processing NMR Data 
 * Copyright (C) 2004-2017 One Moon Scientific, Inc., Westfield, N.J., USA
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package org.nmrfx.processor.datasets.peaks;

import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

/**
 *
 * @author Bruce Johnson
 */
public class ResonanceFactory {

    Map<Long, Resonance> map = new HashMap<>();
    private long lastID = -1;

    public void init() {

    }

    public void clean() {
        Map<Long, Resonance> resonancesNew = new TreeMap<>();
        for (Map.Entry<Long, Resonance> entry : map.entrySet()) {
            Resonance resonance = entry.getValue();
            if (((resonance.getPeakDims() != null) && (!resonance.getPeakDims().isEmpty()))) {
                resonancesNew.put(entry.getKey(), resonance);
            }
        }
        map.clear();
        map = resonancesNew;
    }

    public Resonance build() {
        while (map.get(lastID++)!=null);
        Resonance resonance = new SimpleResonance(lastID);
        map.put(lastID, resonance);
        return resonance;
    }

    public Resonance build(long id) {
        Resonance resonance = get(id);
        if (resonance == null) {
            resonance = new SimpleResonance(id);
            map.put(id,resonance);
        }
        return resonance;
    }

    public Resonance get(long id) {
        return map.get(id);
    }
}
