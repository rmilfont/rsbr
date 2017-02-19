import org.junit.*;
import sbr.*;
import util.*;

public class RegularGraphSegmentationTest {
    @Test
    public void unityGraph() throws Exception {
        Graph noc = RegularGraphBuilder.generateGraph(1, 1);
        testLinearGraph(noc);
    }

    @Test
    public void testLinearGraphs() throws Exception {
        for (int numberOfVertices = 1; numberOfVertices < 10; numberOfVertices++) {
            testHorizontalGraph(numberOfVertices);
            testVerticalGraph(numberOfVertices);
        }
    }

    private void testHorizontalGraph(int numberOfVertices) {
        Graph noc = RegularGraphBuilder.generateGraph(1, numberOfVertices);
        testLinearGraph(noc);
    }

    private void testVerticalGraph(int numberOfVertices) {
        Graph noc = RegularGraphBuilder.generateGraph(numberOfVertices, 1);
        testLinearGraph(noc);
    }

    private void testLinearGraph(Graph noc) {
        int numberOfVertices = noc.getVertices().size();
        SR sr = new SR(noc, new BidimensionalSBRPolicy(noc));
        sr.computeSegments();
        sr.setrestrictions();
        Assert.assertEquals(0, sr.segments().size());
        Assert.assertEquals(0, sr.restrictions().size());
        Assert.assertEquals(numberOfVertices, sr.startVertices().size());
        Assert.assertEquals(numberOfVertices, sr.terminalVertices().size());
    }

    @Test
    public void regular2x2Graph() throws Exception {
        Graph noc = RegularGraphBuilder.generateGraph(2, 2);
        SR sr = new SR(noc, new BidimensionalSBRPolicy(noc));
        sr.computeSegments();
        sr.setrestrictions();
        Assert.assertEquals(1, sr.segments().size());
        Assert.assertEquals(1, sr.restrictions().size());
        Assert.assertEquals(1, sr.startVertices().size());
        Assert.assertEquals(0, sr.terminalVertices().size());
    }

    @Test
    public void regular2x3Graph() throws Exception {
        Graph noc = RegularGraphBuilder.generateGraph(2, 3);
        SR sr = new SR(noc, new BidimensionalSBRPolicy(noc));
        sr.computeSegments();
        sr.setrestrictions();
        Assert.assertEquals(2, sr.segments().size());
        Assert.assertEquals(2, sr.restrictions().size());
        Assert.assertEquals(1, sr.startVertices().size());
        Assert.assertEquals(0, sr.terminalVertices().size());

        SR sr2 = new SR(noc, new BidimensionalSBRPolicy(noc));
        sr2.computeSegments();
        sr2.setrestrictions();

        Assert.assertTrue(sr.restrictions().equals(sr2.restrictions()));
    }

    @Test
    public void regular3x2Graph() throws Exception {
        Graph noc = RegularGraphBuilder.generateGraph(3, 2);
        SR sr = new SR(noc, new BidimensionalSBRPolicy(noc));
        sr.computeSegments();
        sr.setrestrictions();
        Assert.assertEquals(2, sr.segments().size());
        Assert.assertEquals(2, sr.restrictions().size());
        Assert.assertEquals(1, sr.startVertices().size());
        Assert.assertEquals(0, sr.terminalVertices().size());
    }

    @Test
    public void regular3x3Graph() throws Exception {
        Graph noc = RegularGraphBuilder.generateGraph(3, 3);
        SR sr = new SR(noc, new BidimensionalSBRPolicy(noc));
        sr.computeSegments();
        sr.setrestrictions();
        Assert.assertEquals(4, sr.segments().size());
        Assert.assertEquals(4, sr.restrictions().size());
        Assert.assertEquals(1, sr.startVertices().size());
        Assert.assertEquals(0, sr.terminalVertices().size());
    }
}
