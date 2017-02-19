package paths_gen;

import java.util.*;
import util.*;

public class MinimalPathSelector {

    private Graph graph;
    private ArrayList<ArrayList<Path>> paths;
    private LinkWeightTracker linkWeightTracker;

    public MinimalPathSelector(ArrayList<ArrayList<Path>> paths,
                              LinkWeightTracker tracker, Graph graph) {
        this.graph = graph;
        this.paths = paths;
        linkWeightTracker = tracker;
    }

    public ArrayList<ArrayList<Path>> selection() {
        ArrayList<ArrayList<Path>> result = new ArrayList<>();
        for(ArrayList<Path> samePairPaths: paths) {
            result.add(selectPath(samePairPaths));
        }
        return result;
    }

    private ArrayList<Path> selectPath(ArrayList<Path> samePairPaths) {
        Collections.sort(samePairPaths,  new ByLatency(graph));
        Path randomPath = samePairPaths.get(0);
        linkWeightTracker.add(randomPath);
        ArrayList<Path> result = new ArrayList<>();
        result.add(randomPath);
        return result;
    }

    private class ByLatency implements Comparator<Path> {
        Graph graph;

        public ByLatency(Graph graph){
            this.graph = graph;
        }

        @Override
        public int compare(Path p0, Path p1) {
            if(latencyOf(graph, p0) < latencyOf(graph, p1)) return -1;
            if(latencyOf(graph, p0) > latencyOf(graph, p1)) return +1;
            return 0;
        }
    }

    public static Double latencyOf(Graph graph, Path path){
        List<Double> TC_i = new ArrayList<>();
        TC_i.add(1.0);
        for(int i = 0; i < path.size()-1; i++) {
            Edge e = graph.adjunct(path.get(i), path.get(i+1));
            TC_i.add(1/(1-e.weight()));
        }
        TC_i.add(1.0);

        JMMatlab jmMatlab = new JMMatlab();
        return jmMatlab.latencia(32, TC_i);
    }
}