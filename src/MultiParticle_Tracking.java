import ij.*;
import ij.plugin.PlugIn;
import ij.process.*;
import ij.gui.*;
import ij.plugin.filter.ThresholdToSelection;
import ij.plugin.frame.RoiManager;
import ij.measure.ResultsTable;
import ij.plugin.filter.ParticleAnalyzer;
import ij.measure.Calibration;
import ij.measure.Measurements;

import java.awt.*;
import java.awt.geom.PathIterator;
import java.awt.geom.Point2D;
import java.util.*;
import javax.swing.*;
import javax.swing.border.*;

/**
 * Multi-Particle Tracking Plugin
 * Function: Track particle detachment time within probe region and generate pseudo-color map
 * Version: 7.1 (Fixed variable error)
 */
public class MultiParticle_Tracking implements PlugIn {
    
    /* ==================== Parameter Settings ==================== */
    private double pixelSize = 1.0;           // Pixel size (microns/pixel)
    private String unit = "μm";                // Unit
    private double frameInterval = 0.0;        // Frame interval (milliseconds)
    private Roi userROI = null;                // User drawn probe region (pentagon) - not added to ROI Manager
    private Polygon probePolygon = null;       // Probe polygon
    private Point2D.Double probeCentroid = null; // Probe centroid
    private double[] probeToEdgeDistances = new double[5]; // Shortest distances from probe centroid to five edges (pixels)
    
    /* ==================== Analysis Parameters ==================== */
    private double minArea = 10.0;              // Minimum area (μm²)
    private double maxArea = 1000.0;             // Maximum area (μm²)
    private double minCircularity = 0.5;         // Minimum circularity
    private double maxCircularity = 1.0;         // Maximum circularity
    private double distanceDeviation = 10.0;     // Cell distance deviation threshold (default 10%)
    
    /* ==================== Image Related ==================== */
    private ImagePlus currentImp = null;         // Current image
    private ImagePlus previewWindow = null;       // Preview window
    private ImagePlus pseudoColorWindow = null;   // Pseudo-color window
    private String baseResultsTitle = "Probe Region Cell Detachment Analysis Results"; // Base results table title
    
    /* ==================== Result Storage ==================== */
    private Map<Integer, ParticleData> particles = new HashMap<Integer, ParticleData>();
    private double maxDropTime = 0;               // Maximum detachment time
    
    /* ==================== Pseudo-color Settings ==================== */
    private String currentLutType = "Jet";        // Current LUT type
    private double currentMinThreshold = 0;        // Current minimum threshold
    private double currentMaxThreshold = 0;        // Current maximum threshold
    private double detachmentConstant = 10.78;     // Detachment constant C (default 10.78)
    private JFrame pseudoFrame;                    // Pseudo-color settings window
    private JPanel lutPreviewPanel;                 // LUT preview panel
    
    /* ==================== UI Components ==================== */
    private JFrame mainFrame;
    private JTextField fileField;                // File name display
    private JTextField scaleField;                // Scale display (μm/pixel)
    private JTextField frameTimeField;            // Frame interval (milliseconds)
    private JLabel roiStatusLabel;                // Probe status
    private JTextField minAreaField;              // Minimum area input (μm²)
    private JTextField maxAreaField;              // Maximum area input (μm²)
    private JTextField minCircularityField;        // Minimum circularity input
    private JTextField maxCircularityField;        // Maximum circularity input
    private JTextField distanceDeviationField;     // Cell distance deviation percentage
    private JButton previewButton;                 // Preview button
    private JButton analyzeButton;                 // Analyze button
    private JTextArea logArea;                     // Log area
    private JPanel mainPanel;                       // Main panel
    private JDialog summaryDialog;                  // Summary dialog
    
    // 新策略：使用静态计数器来生成唯一的表格标题
    private static int analysisCounter = 0;
    
    /* ==================== Inner Class: Particle Data Storage ==================== */
    class ParticleData {
        int id;                          // Particle ID
        double initialArea;               // Initial area (μm²)
        double initialPerimeter;           // Initial perimeter (μm)
        double initialCircularity;         // Initial circularity
        double distanceDeviation;          // Distance deviation from probe (%)
        double dropTime = -1;              // Detachment time (seconds), -1 means not detached
        int dropFrame = -1;                // Detachment frame number
        Roi initialRoi;                    // Initial contour
        Point2D.Double centroid;            // Cell centroid (pixel coordinates)
        
        ParticleData(int id, double area, double perimeter, double circularity, 
                    double deviation, Point2D.Double centroid, Roi roi) {
            this.id = id;
            this.initialArea = area;
            this.initialPerimeter = perimeter;
            this.initialCircularity = circularity;
            this.distanceDeviation = deviation;
            this.centroid = centroid;
            this.initialRoi = roi;
        }
    }
    
    /* ==================== Geometry Calculation Methods ==================== */
    
    /**
     * Calculate shortest distance from point to line segment (considering endpoints)
     */
    public static double pointToSegmentDistance(double px, double py,
                                               double x1, double y1,
                                               double x2, double y2) {
        double dx = x2 - x1;
        double dy = y2 - y1;
        
        // If segment degenerates to a point
        if (dx == 0 && dy == 0) {
            return Math.hypot(px - x1, py - y1);
        }
        
        // Calculate projection parameter t
        double t = ((px - x1) * dx + (py - y1) * dy) / (dx * dx + dy * dy);
        
        if (t < 0) {
            // Projection outside segment start
            return Math.hypot(px - x1, py - y1);
        } else if (t > 1) {
            // Projection outside segment end
            return Math.hypot(px - x2, py - y2);
        } else {
            // Projection on segment
            double projX = x1 + t * dx;
            double projY = y1 + t * dy;
            return Math.hypot(px - projX, py - projY);
        }
    }
    
    /**
     * Check if two line segments intersect
     */
    public static boolean segmentsIntersect(double x1, double y1, double x2, double y2,
                                           double x3, double y3, double x4, double y4) {
        double d1 = direction(x3, y3, x4, y4, x1, y1);
        double d2 = direction(x3, y3, x4, y4, x2, y2);
        double d3 = direction(x1, y1, x2, y2, x3, y3);
        double d4 = direction(x1, y1, x2, y2, x4, y4);
        
        if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
            ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) {
            return true;
        }
        
        if (d1 == 0 && onSegment(x3, y3, x4, y4, x1, y1)) return true;
        if (d2 == 0 && onSegment(x3, y3, x4, y4, x2, y2)) return true;
        if (d3 == 0 && onSegment(x1, y1, x2, y2, x3, y3)) return true;
        if (d4 == 0 && onSegment(x1, y1, x2, y2, x4, y4)) return true;
        
        return false;
    }
    
    private static double direction(double x1, double y1, double x2, double y2, double x3, double y3) {
        return (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
    }
    
    private static boolean onSegment(double x1, double y1, double x2, double y2, double x, double y) {
        return Math.min(x1, x2) <= x && x <= Math.max(x1, x2) &&
               Math.min(y1, y2) <= y && y <= Math.max(y1, y2);
    }
    
    /**
     * Check if cell contour intersects any edge of the pentagon
     */
    public static boolean cellIntersectsProbe(Roi cellRoi, Polygon probePoly) {
        Polygon cellPoly = cellRoi.getPolygon();
        int nCell = cellPoly.npoints;
        int nProbe = probePoly.npoints;
        
        for (int i = 0; i < nCell; i++) {
            double cx1 = cellPoly.xpoints[i];
            double cy1 = cellPoly.ypoints[i];
            double cx2 = cellPoly.xpoints[(i + 1) % nCell];
            double cy2 = cellPoly.ypoints[(i + 1) % nCell];
            
            for (int j = 0; j < nProbe; j++) {
                double px1 = probePoly.xpoints[j];
                double py1 = probePoly.ypoints[j];
                double px2 = probePoly.xpoints[(j + 1) % nProbe];
                double py2 = probePoly.ypoints[(j + 1) % nProbe];
                
                if (segmentsIntersect(cx1, cy1, cx2, cy2, px1, py1, px2, py2)) {
                    return true;
                }
            }
        }
        return false;
    }
    
    /**
     * Calculate polygon centroid
     */
    public static Point2D.Double calculatePolygonCentroid(Polygon poly) {
        double area = 0;
        double cx = 0;
        double cy = 0;
        
        int n = poly.npoints;
        for (int i = 0; i < n; i++) {
            int x1 = poly.xpoints[i];
            int y1 = poly.ypoints[i];
            int x2 = poly.xpoints[(i + 1) % n];
            int y2 = poly.ypoints[(i + 1) % n];
            
            double cross = x1 * y2 - x2 * y1;
            area += cross;
            cx += (x1 + x2) * cross;
            cy += (y1 + y2) * cross;
        }
        
        area /= 2.0;
        cx /= (6.0 * area);
        cy /= (6.0 * area);
        
        return new Point2D.Double(cx, cy);
    }
    
    /**
     * Calculate distances from probe centroid to five edges
     */
    private double[] calculateProbeToEdgeDistances(Polygon probePoly, Point2D.Double centroid) {
        double[] distances = new double[5];
        
        for (int i = 0; i < probePoly.npoints; i++) {
            int x1 = probePoly.xpoints[i];
            int y1 = probePoly.ypoints[i];
            int x2 = probePoly.xpoints[(i + 1) % probePoly.npoints];
            int y2 = probePoly.ypoints[(i + 1) % probePoly.npoints];
            
            distances[i] = pointToSegmentDistance(centroid.x, centroid.y, x1, y1, x2, y2);
        }
        
        return distances;
    }
    
    /**
     * Calculate minimum distance deviation from cell to probe
     */
    private double calculateDeviation(Roi cellRoi, Point2D.Double cellCentroid, 
                                      Polygon probePoly, Point2D.Double probeCentroid,
                                      double[] probeToEdgeDistances) {
        
        if (cellIntersectsProbe(cellRoi, probePoly)) {
            return 0.0;
        }
        
        double[] cellToEdgeDistances = new double[probePoly.npoints];
        
        for (int i = 0; i < probePoly.npoints; i++) {
            int x1 = probePoly.xpoints[i];
            int y1 = probePoly.ypoints[i];
            int x2 = probePoly.xpoints[(i + 1) % probePoly.npoints];
            int y2 = probePoly.ypoints[(i + 1) % probePoly.npoints];
            
            cellToEdgeDistances[i] = pointToSegmentDistance(cellCentroid.x, cellCentroid.y, x1, y1, x2, y2);
        }
        
        double minCellDist = Double.MAX_VALUE;
        int minEdgeIndex = -1;
        
        for (int i = 0; i < cellToEdgeDistances.length; i++) {
            if (cellToEdgeDistances[i] < minCellDist) {
                minCellDist = cellToEdgeDistances[i];
                minEdgeIndex = i;
            }
        }
        
        double probeDist = probeToEdgeDistances[minEdgeIndex];
        
        if (probeDist == 0) return 100.0;
        return (minCellDist / probeDist) * 100.0;
    }
    
    /* ==================== Plugin Main Method ==================== */
    @Override
    public void run(String arg) {
        createIntegratedUI();
    }
    
    /* ==================== UI Creation ==================== */
    private void createIntegratedUI() {
        mainFrame = new JFrame("Multi-Particle Tracking Analysis");
        mainFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        mainFrame.addWindowListener(new java.awt.event.WindowAdapter() {
            @Override
            public void windowClosing(java.awt.event.WindowEvent windowEvent) {
                closeAllWindows();
            }
        });
        
        mainPanel = new JPanel();
        mainPanel.setLayout(new BoxLayout(mainPanel, BoxLayout.Y_AXIS));
        mainPanel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
        
        /* ---------- 1. File Information Panel ---------- */
        JPanel filePanel = createLabeledPanel("File Information");
        fileField = new JTextField("No file opened", 25);
        fileField.setEditable(false);
        JButton openButton = new JButton("Open Video/Sequence");
        openButton.addActionListener(e -> openVideoFile());
        
        JPanel fileControlPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 5, 0));
        fileControlPanel.add(new JLabel("Current file:"));
        fileControlPanel.add(fileField);
        fileControlPanel.add(openButton);
        filePanel.add(fileControlPanel);
        
        /* ---------- 2. Scale Setting Panel ---------- */
        JPanel scalePanel = createLabeledPanel("Scale Setting");
        scalePanel.setLayout(new FlowLayout(FlowLayout.LEFT, 5, 0));
        
        JLabel currentScaleLabel = new JLabel("Scale:");
        scaleField = new JTextField("Not set", 8);
        scaleField.setEditable(false);
        scaleField.setBackground(Color.LIGHT_GRAY);
        
        JLabel unitLabel = new JLabel("μm/pixel");
        
        JButton setScaleButton = new JButton("Set");
        setScaleButton.addActionListener(e -> setupScale());
        
        scalePanel.add(currentScaleLabel);
        scalePanel.add(scaleField);
        scalePanel.add(unitLabel);
        scalePanel.add(setScaleButton);
        
        /* ---------- 3. Time Setting Panel ---------- */
        JPanel timePanel = createLabeledPanel("Time Setting");
        JLabel timeLabel = new JLabel("Frame interval (ms):");
        frameTimeField = new JTextField("200", 10);
        
        JPanel timeControlPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 5, 0));
        timeControlPanel.add(timeLabel);
        timeControlPanel.add(frameTimeField);
        timePanel.add(timeControlPanel);
        
        /* ---------- 4. Probe Region Setting Panel ---------- */
        JPanel roiPanel = createLabeledPanel("Probe Region Setting");
        roiStatusLabel = new JLabel("Current probe region: Not set");
        
        JPanel roiButtonPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 5, 0));
        JButton drawRoiButton = new JButton("Draw Pentagon Probe");
        drawRoiButton.addActionListener(e -> drawPentagonROI());
        
        JButton loadRoiButton = new JButton("Load Existing ROI");
        loadRoiButton.addActionListener(e -> loadExistingROI());
        
        roiButtonPanel.add(drawRoiButton);
        roiButtonPanel.add(loadRoiButton);
        
        roiPanel.add(roiStatusLabel, BorderLayout.NORTH);
        roiPanel.add(roiButtonPanel, BorderLayout.CENTER);
        
        /* ---------- 5. Cell Detection Parameters Panel ---------- */
        JPanel detectionPanel = createLabeledPanel("Cell Detection Parameters");
        detectionPanel.setLayout(new GridLayout(5, 2, 5, 5));
        
        minAreaField = new JTextField("10.0", 8);
        maxAreaField = new JTextField("1000.0", 8);
        minCircularityField = new JTextField("0.5", 8);
        maxCircularityField = new JTextField("1.0", 8);
        distanceDeviationField = new JTextField("10.0", 8);
        
        detectionPanel.add(new JLabel("Min Area (μm²):"));
        detectionPanel.add(minAreaField);
        detectionPanel.add(new JLabel("Max Area (μm²):"));
        detectionPanel.add(maxAreaField);
        detectionPanel.add(new JLabel("Min Circularity (0-1):"));
        detectionPanel.add(minCircularityField);
        detectionPanel.add(new JLabel("Max Circularity (0-1):"));
        detectionPanel.add(maxCircularityField);
        detectionPanel.add(new JLabel("Cell Distance Deviation (%):"));
        detectionPanel.add(distanceDeviationField);
        
        /* ---------- 6. Action Buttons Panel ---------- */
        JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.CENTER, 10, 5));
        previewButton = new JButton("Preview First Frame");
        previewButton.addActionListener(e -> previewFirstFrame());
        
        analyzeButton = new JButton("Start Tracking Analysis");
        analyzeButton.addActionListener(e -> startAnalysis());
        
        buttonPanel.add(previewButton);
        buttonPanel.add(analyzeButton);
        
        /* ---------- 7. Log Panel ---------- */
        JPanel logPanel = createLabeledPanel("Log");
        logArea = new JTextArea(6, 50);
        logArea.setEditable(false);
        JScrollPane scrollPane = new JScrollPane(logArea);
        logPanel.add(scrollPane);
        
        /* ---------- Assemble all panels ---------- */
        mainPanel.add(filePanel);
        mainPanel.add(Box.createVerticalStrut(8));
        mainPanel.add(scalePanel);
        mainPanel.add(Box.createVerticalStrut(8));
        mainPanel.add(timePanel);
        mainPanel.add(Box.createVerticalStrut(8));
        mainPanel.add(roiPanel);
        mainPanel.add(Box.createVerticalStrut(8));
        mainPanel.add(detectionPanel);
        mainPanel.add(Box.createVerticalStrut(8));
        mainPanel.add(buttonPanel);
        mainPanel.add(Box.createVerticalStrut(8));
        mainPanel.add(logPanel);
        
        mainFrame.add(mainPanel, BorderLayout.CENTER);
        mainFrame.pack();
        mainFrame.setLocationRelativeTo(null);
        mainFrame.setVisible(true);
        
        log("Multi-Particle Tracking Plugin started");
    }
    
    private JPanel createLabeledPanel(String title) {
        JPanel panel = new JPanel();
        panel.setLayout(new BorderLayout());
        panel.setBorder(BorderFactory.createTitledBorder(
            BorderFactory.createEtchedBorder(), title));
        return panel;
    }
    
    private void log(String message) {
        String timestamp = new java.text.SimpleDateFormat("HH:mm:ss").format(new java.util.Date());
        logArea.append("[" + timestamp + "] " + message + "\n");
        logArea.setCaretPosition(logArea.getDocument().getLength());
    }
    
    /* ==================== Close all related windows ==================== */
    private void closeAllWindows() {
        // Close preview window
        if (previewWindow != null) {
            previewWindow.close();
            previewWindow = null;
        }
        
        // Close pseudo-color window
        if (pseudoColorWindow != null) {
            pseudoColorWindow.close();
            pseudoColorWindow = null;
        }
        
        // Close pseudo-color settings window
        if (pseudoFrame != null) {
            pseudoFrame.dispose();
            pseudoFrame = null;
        }
        
        // Close summary dialog
        if (summaryDialog != null) {
            summaryDialog.dispose();
            summaryDialog = null;
        }
        
        // Close all result windows - 使用 baseResultsTitle 替代 currentResultsTitle
        Frame[] frames = Frame.getFrames();
        for (Frame frame : frames) {
            String title = frame.getTitle();
            if (title != null && title.contains(baseResultsTitle)) {
                frame.dispose();
            }
        }
        
        // Close all ImageJ windows (except main interface)
        int[] windowIDs = WindowManager.getIDList();
        if (windowIDs != null) {
            for (int id : windowIDs) {
                ImagePlus imp = WindowManager.getImage(id);
                if (imp != null && imp != currentImp) {
                    imp.close();
                }
            }
        }
        
        // Close current image window
        if (currentImp != null) {
            currentImp.close();
            currentImp = null;
        }
        
        // Clear ROI Manager
        RoiManager roiManager = RoiManager.getInstance();
        if (roiManager != null) {
            roiManager.reset();
        }
        
        // Reset state
        particles.clear();
        userROI = null;
        probePolygon = null;
        probeCentroid = null;
        probeToEdgeDistances = null;
        maxDropTime = 0;
        
        // Reset UI
        fileField.setText("No file opened");
        scaleField.setText("Not set");
        roiStatusLabel.setText("Current probe region: Not set");
        pixelSize = 1.0;
    }
    
    /* ==================== Function 1: Open video file ==================== */
    private void openVideoFile() {
        SwingWorker<ImagePlus, String> worker = new SwingWorker<ImagePlus, String>() {
            @Override
            protected ImagePlus doInBackground() throws Exception {
                publish("Opening file dialog...");
                return IJ.openImage();
            }
            
            @Override
            protected void process(java.util.List<String> chunks) {
                for (String msg : chunks) {
                    log(msg);
                }
            }
            
            @Override
            protected void done() {
                try {
                    ImagePlus imp = get();
                    if (imp != null) {
                        // Close all related windows
                        closeAllWindows();
                        
                        // Reset analysis counter for new file
                        analysisCounter = 0;
                        
                        // Set new image
                        currentImp = imp;
                        new Thread(() -> {
                            currentImp.show();
                        }).start();
                        
                        fileField.setText(currentImp.getTitle());
                        
                        log("File opened: " + currentImp.getTitle());
                        log("Image size: " + currentImp.getWidth() + "x" + currentImp.getHeight());
                        log("Total frames: " + currentImp.getStackSize());
                    } else {
                        log("No file selected");
                    }
                } catch (Exception e) {
                    log("Failed to open file: " + e.getMessage());
                    e.printStackTrace();
                }
            }
        };
        worker.execute();
    }
    
    /* ==================== Function 2: Set scale ==================== */
    private void setupScale() {
        if (currentImp == null) {
            IJ.showMessage("Please open an image first");
            return;
        }
        
        log("Please draw a line of known length using the line tool, then click Confirm");
        
        Toolbar.getInstance().setTool(Toolbar.LINE);
        WindowManager.setCurrentWindow(currentImp.getWindow());
        
        JDialog hintDialog = new JDialog(mainFrame, "Set Scale", false);
        hintDialog.setLayout(new BorderLayout());
        hintDialog.setSize(300, 180);
        hintDialog.setLocationRelativeTo(mainFrame);
        
        JPanel hintPanel = new JPanel();
        hintPanel.setLayout(new BoxLayout(hintPanel, BoxLayout.Y_AXIS));
        hintPanel.setBorder(BorderFactory.createEmptyBorder(10, 10, 20, 10));
        
        JLabel hintLabel = new JLabel("<html><center>Please draw a line of known length using the line tool<br>Click 'Confirm' when done</center></html>");
        hintLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
        hintLabel.setFont(new Font("SansSerif", Font.PLAIN, 12));
        
        JButton okButton = new JButton("Confirm");
        okButton.setAlignmentX(Component.CENTER_ALIGNMENT);
        okButton.addActionListener(ev -> {
            hintDialog.dispose();
            
            Roi roi = currentImp.getRoi();
            if (roi == null || !(roi instanceof Line)) {
                IJ.error("No line detected, please draw a line first");
                return;
            }
            
            double pixelLength = ((Line) roi).getLength();
            
            // Create custom length input dialog with consistent style
            JDialog lengthDialog = new JDialog(mainFrame, "Set Scale", true);
            lengthDialog.setLayout(new BorderLayout());
            lengthDialog.setSize(300, 180);
            lengthDialog.setLocationRelativeTo(mainFrame);
            
            JPanel lengthPanel = new JPanel();
            lengthPanel.setLayout(new BoxLayout(lengthPanel, BoxLayout.Y_AXIS));
            lengthPanel.setBorder(BorderFactory.createEmptyBorder(10, 10, 20, 10));
            
            JLabel lengthLabel = new JLabel("<html><center>Please enter the actual length of the drawn line<br>Unit: μm</center></html>");
            lengthLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
            lengthLabel.setFont(new Font("SansSerif", Font.PLAIN, 12));
            
            JTextField lengthField = new JTextField("100.0", 10);
            lengthField.setMaximumSize(new Dimension(200, 25));
            lengthField.setAlignmentX(Component.CENTER_ALIGNMENT);
            
            JPanel lengthButtonPanel = new JPanel(new FlowLayout(FlowLayout.CENTER));
            JButton lengthOkButton = new JButton("Confirm");
            JButton lengthCancelButton = new JButton("Cancel");
            
            lengthOkButton.addActionListener(ev2 -> {
                try {
                    double realLength = Double.parseDouble(lengthField.getText());
                    pixelSize = realLength / pixelLength;
                    scaleField.setText(String.format("%.4f", pixelSize));
                    
                    Calibration cal = currentImp.getCalibration();
                    cal.pixelWidth = pixelSize;
                    cal.pixelHeight = pixelSize;
                    cal.setUnit("μm");
                    
                    log(String.format("Scale set: %.4f μm/pixel (%.1f pixels = %.1f μm)", 
                        pixelSize, pixelLength, realLength));
                    
                    lengthDialog.dispose();
                } catch (NumberFormatException ex) {
                    IJ.error("Please enter a valid number");
                }
            });
            
            lengthCancelButton.addActionListener(ev2 -> {
                lengthDialog.dispose();
                log("Scale setting canceled");
            });
            
            lengthButtonPanel.add(lengthOkButton);
            lengthButtonPanel.add(lengthCancelButton);
            
            lengthPanel.add(lengthLabel);
            lengthPanel.add(Box.createVerticalStrut(10));
            lengthPanel.add(lengthField);
            lengthPanel.add(Box.createVerticalStrut(20));
            lengthPanel.add(lengthButtonPanel);
            
            lengthDialog.add(lengthPanel);
            lengthDialog.setVisible(true);
        });
        
        JButton cancelButton = new JButton("Cancel");
        cancelButton.setAlignmentX(Component.CENTER_ALIGNMENT);
        cancelButton.addActionListener(ev -> {
            hintDialog.dispose();
            log("Scale setting canceled");
        });
        
        JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.CENTER));
        buttonPanel.add(okButton);
        buttonPanel.add(cancelButton);
        
        hintPanel.add(hintLabel);
        hintPanel.add(Box.createVerticalStrut(20));
        hintPanel.add(buttonPanel);
        
        hintDialog.add(hintPanel);
        hintDialog.setVisible(true);
    }
    
    /* ==================== Function 3: Draw pentagon probe (not added to ROI Manager) ==================== */
    private void drawPentagonROI() {
        if (currentImp == null) {
            IJ.showMessage("Please open an image first");
            return;
        }
        
        log("Please draw a pentagon probe region using the polygon tool, then click Confirm");
        
        Toolbar.getInstance().setTool(Toolbar.POLYGON);
        WindowManager.setCurrentWindow(currentImp.getWindow());
        
        JDialog hintDialog = new JDialog(mainFrame, "Draw Pentagon Probe", false);
        hintDialog.setLayout(new BorderLayout());
        hintDialog.setSize(350, 230);
        hintDialog.setLocationRelativeTo(mainFrame);
        
        JPanel hintPanel = new JPanel();
        hintPanel.setLayout(new BoxLayout(hintPanel, BoxLayout.Y_AXIS));
        hintPanel.setBorder(BorderFactory.createEmptyBorder(10, 10, 20, 10));
        
        JLabel hintLabel = new JLabel("<html><center>Please draw a pentagon probe region using the polygon tool<br>Steps:<br>1. Click 5 points to form a pentagon<br>2. Press Enter to complete<br>3. Click the 'Confirm' button below</center></html>");
        hintLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
        hintLabel.setFont(new Font("SansSerif", Font.PLAIN, 12));
        
        JButton okButton = new JButton("Confirm");
        okButton.setAlignmentX(Component.CENTER_ALIGNMENT);
        okButton.addActionListener(ev -> {
            hintDialog.dispose();
            
            Roi roi = currentImp.getRoi();
            if (roi == null) {
                IJ.error("No probe region detected, please draw a pentagon first");
                return;
            }
            
            if (roi.getType() == Roi.POLYGON) {
                Polygon poly = roi.getPolygon();
                if (poly.npoints == 5) {
                    // Directly save probe, not added to ROI Manager
                    userROI = roi;
                    probePolygon = poly;
                    probeCentroid = calculatePolygonCentroid(poly);
                    probeToEdgeDistances = calculateProbeToEdgeDistances(poly, probeCentroid);
                    
                    roiStatusLabel.setText("Current probe region: Set (pentagon)");
                    log("Pentagon probe region set");
                    log("Probe centroid: (" + String.format("%.1f", probeCentroid.x) + ", " + 
                        String.format("%.1f", probeCentroid.y) + ") pixels");
                    
                    for (int i = 0; i < 5; i++) {
                        log(String.format("Edge %d distance: %.1f pixels (%.2f μm)", 
                            i+1, probeToEdgeDistances[i], probeToEdgeDistances[i] * pixelSize));
                    }
                } else {
                    IJ.error("Probe region must be a pentagon, you drew a " + poly.npoints + "-gon");
                }
            } else {
                IJ.error("Please use the polygon tool to draw a pentagon region");
            }
        });
        
        JButton cancelButton = new JButton("Cancel");
        cancelButton.setAlignmentX(Component.CENTER_ALIGNMENT);
        cancelButton.addActionListener(ev -> {
            hintDialog.dispose();
            log("Probe drawing canceled");
        });
        
        JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.CENTER));
        buttonPanel.add(okButton);
        buttonPanel.add(cancelButton);
        
        hintPanel.add(hintLabel);
        hintPanel.add(Box.createVerticalStrut(20));
        hintPanel.add(buttonPanel);
        
        hintDialog.add(hintPanel);
        hintDialog.setVisible(true);
    }
    
    /* ==================== Function 4: Load existing ROI ==================== */
    private void loadExistingROI() {
        if (currentImp == null) {
            IJ.showMessage("Please open an image first");
            return;
        }
        
        // Get ROI Manager - use final array to solve anonymous inner class access issue
        final RoiManager[] roiManagerHolder = new RoiManager[1];
        roiManagerHolder[0] = RoiManager.getInstance();
        if (roiManagerHolder[0] == null) {
            roiManagerHolder[0] = new RoiManager();
        }
        
        final RoiManager roiManager = roiManagerHolder[0];
        roiManager.setVisible(true);
        roiManager.runCommand("Show All");
        
        log("Please select a pentagon ROI in the ROI Manager as probe region");
        
        // Create non-modal hint dialog
        JDialog hintDialog = new JDialog(mainFrame, "Select ROI", false);
        hintDialog.setLayout(new BorderLayout());
        hintDialog.setSize(350, 180);
        hintDialog.setLocationRelativeTo(mainFrame);
        
        JPanel hintPanel = new JPanel();
        hintPanel.setLayout(new BoxLayout(hintPanel, BoxLayout.Y_AXIS));
        hintPanel.setBorder(BorderFactory.createEmptyBorder(10, 10, 20, 10));
        
        JLabel hintLabel = new JLabel("<html><center>Please select a pentagon ROI in the ROI Manager<br>Click 'Confirm' after selection</center></html>");
        hintLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
        hintLabel.setFont(new Font("SansSerif", Font.PLAIN, 12));
        
        JButton okButton = new JButton("Confirm");
        okButton.setAlignmentX(Component.CENTER_ALIGNMENT);
        okButton.addActionListener(ev -> {
            hintDialog.dispose();
            
            // Get selected ROI from ROI Manager
            Roi[] rois = roiManager.getSelectedRoisAsArray();
            if (rois.length > 0) {
                Roi selectedRoi = rois[0];
                // Check if it's a pentagon
                if (selectedRoi.getType() == Roi.POLYGON) {
                    Polygon poly = selectedRoi.getPolygon();
                    if (poly.npoints == 5) {
                        userROI = selectedRoi;
                        probePolygon = poly;
                        probeCentroid = calculatePolygonCentroid(poly);
                        probeToEdgeDistances = calculateProbeToEdgeDistances(poly, probeCentroid);
                        
                        roiStatusLabel.setText("Current probe region: Loaded (pentagon)");
                        log("Pentagon probe region loaded");
                        log("Probe centroid: (" + String.format("%.1f", probeCentroid.x) + ", " + 
                            String.format("%.1f", probeCentroid.y) + ") pixels");
                        
                        for (int i = 0; i < 5; i++) {
                            log(String.format("Edge %d distance: %.1f pixels (%.2f μm)", 
                                i+1, probeToEdgeDistances[i], probeToEdgeDistances[i] * pixelSize));
                        }
                    } else {
                        IJ.error("Selected ROI must be a pentagon");
                    }
                } else {
                    IJ.error("Please select a polygon type ROI");
                }
            } else {
                IJ.error("Please select an ROI in the ROI Manager first");
            }
        });
        
        JButton cancelButton = new JButton("Cancel");
        cancelButton.setAlignmentX(Component.CENTER_ALIGNMENT);
        cancelButton.addActionListener(ev -> {
            hintDialog.dispose();
            log("ROI loading canceled");
        });
        
        JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.CENTER));
        buttonPanel.add(okButton);
        buttonPanel.add(cancelButton);
        
        hintPanel.add(hintLabel);
        hintPanel.add(Box.createVerticalStrut(20));
        hintPanel.add(buttonPanel);
        
        hintDialog.add(hintPanel);
        hintDialog.setVisible(true);
    }
    
    /* ==================== Function 5: Preview first frame detection ==================== */
    private void previewFirstFrame() {
        if (!validateParameters()) {
            return;
        }
        
        if (previewWindow != null) {
            previewWindow.close();
        }
        
        log("Previewing first frame cell detection...");
        log(String.format("Current scale: %.4f μm/pixel", pixelSize));
        log(String.format("Area threshold: %.1f - %.1f μm² (pixels: %.1f - %.1f px²)", 
            minArea, maxArea, 
            minArea / (pixelSize * pixelSize), 
            maxArea / (pixelSize * pixelSize)));
        
        SwingWorker<Void, String> worker = new SwingWorker<Void, String>() {
            @Override
            protected Void doInBackground() throws Exception {
                // Clear previous particle data
                particles.clear();
                
                // Clear ROI Manager
                RoiManager roiManager = RoiManager.getInstance();
                if (roiManager != null) {
                    roiManager.reset();
                }
                
                // Execute preview (detect first frame cells)
                detectFirstFrameCells();
                
                publish("First frame detection completed, found " + particles.size() + " eligible cells");
                
                return null;
            }
            
            @Override
            protected void process(java.util.List<String> chunks) {
                for (String msg : chunks) {
                    log(msg);
                }
            }
            
            @Override
            protected void done() {
                try {
                    get();
                } catch (Exception e) {
                    log("Preview failed: " + e.getMessage());
                    e.printStackTrace();
                }
            }
        };
        worker.execute();
    }
    
    /* ==================== Function 6: Start tracking analysis (auto-show pseudo-color) ==================== */
    private void startAnalysis() {
        if (!validateParameters()) {
            return;
        }
        
        // Check if preview has been done
        if (particles.isEmpty()) {
            IJ.showMessage("Please run preview first frame detection first");
            return;
        }
        
        log("Starting tracking analysis...");
        
        SwingWorker<Void, String> worker = new SwingWorker<Void, String>() {
            @Override
            protected Void doInBackground() throws Exception {
                // Directly use cells detected in preview
                // No need to re-detect first frame
                
                publish("Tracking " + particles.size() + " cells...");
                
                // Execute tracking
                trackAllFrames();
                
                // Calculate maximum detachment time
                calculateMaxDropTime();
                
                // Display results
                SwingUtilities.invokeAndWait(new Runnable() {
                    @Override
                    public void run() {
                        // 新策略：每次使用不同的标题
                        analysisCounter++;
                        displayResults();
                        
                        // Then show pseudo-color settings window
                        showPseudoColorFrame();
                        showSummary();
                    }
                });
                
                publish("Analysis completed!");
                
                return null;
            }
            
            @Override
            protected void process(java.util.List<String> chunks) {
                for (String msg : chunks) {
                    log(msg);
                    previewButton.setEnabled(false);
                    analyzeButton.setEnabled(false);
                }
            }
            
            @Override
            protected void done() {
                try {
                    get();
                    previewButton.setEnabled(true);
                    analyzeButton.setEnabled(true);
                } catch (Exception e) {
                    log("Error during analysis: " + e.getMessage());
                    previewButton.setEnabled(true);
                    analyzeButton.setEnabled(true);
                    e.printStackTrace();
                }
            }
        };
        worker.execute();
    }
    
    private boolean validateParameters() {
        if (currentImp == null) {
            IJ.showMessage("Please open an image first");
            return false;
        }
        
        if (scaleField.getText().equals("Not set")) {
            IJ.showMessage("Please set scale first");
            return false;
        }
        
        try {
            frameInterval = Double.parseDouble(frameTimeField.getText());
            minArea = Double.parseDouble(minAreaField.getText());
            maxArea = Double.parseDouble(maxAreaField.getText());
            minCircularity = Double.parseDouble(minCircularityField.getText());
            maxCircularity = Double.parseDouble(maxCircularityField.getText());
            distanceDeviation = Double.parseDouble(distanceDeviationField.getText());
            
            if (frameInterval <= 0) {
                IJ.showMessage("Frame interval must be greater than 0");
                return false;
            }
            
            if (minArea <= 0 || maxArea <= minArea) {
                IJ.showMessage("Invalid area parameters");
                return false;
            }
            
            if (minCircularity < 0 || minCircularity > 1 || 
                maxCircularity < 0 || maxCircularity > 1 || 
                minCircularity >= maxCircularity) {
                IJ.showMessage("Invalid circularity parameters, must be between 0-1 and min < max");
                return false;
            }
            
            if (distanceDeviation < 0 || distanceDeviation > 10000) {
                IJ.showMessage("Cell distance deviation must be between 0-10000%");
                return false;
            }
            
            if (userROI == null) {
                IJ.showMessage("Please set probe region first");
                return false;
            }
            
            return true;
            
        } catch (NumberFormatException e) {
            IJ.showMessage("Please enter valid numbers");
            return false;
        }
    }
    
    /* ==================== Cell segmentation method ==================== */
    private ImageProcessor segmentCells(ImageProcessor ip) {
        ImageProcessor processed = ip.duplicate();
        
        // 1. Convert to 8-bit
        if (!(processed instanceof ByteProcessor)) {
            processed = processed.convertToByte(true);
        }
        
        // 2. Gaussian blur to reduce noise (sigma=1)
        processed.blurGaussian(1.0);
        
        // 3. Enhance contrast
        processed.resetMinAndMax();
        
        // 4. Set auto threshold (Default dark)
        processed.setAutoThreshold("Default dark");
        processed.threshold(processed.getAutoThreshold());
        
        // 5. Morphological operations to improve segmentation
        processed.erode();
        processed.dilate();
        
        // 6. Create binary mask
        ImageProcessor mask = processed.createMask();
        
        return mask;
    }
    
    private void detectFirstFrameCells() {
        ImageProcessor ip = currentImp.getProcessor().duplicate();
        if (!(ip instanceof ByteProcessor)) {
            ip = ip.convertToByte(true);
        }
        
        previewWindow = new ImagePlus("First Frame Cell Detection", ip.duplicate());
        
        // Use reference algorithm segmentation method
        ImageProcessor segmented = segmentCells(ip);
        
        // Get ROI Manager
        RoiManager roiManager = RoiManager.getInstance();
        if (roiManager == null) {
            roiManager = new RoiManager();
        }
        
        // Clear ROI Manager
        roiManager.reset();
        
        // Calculate pixel unit area thresholds
        double minAreaPixels = minArea / (pixelSize * pixelSize);
        double maxAreaPixels = maxArea / (pixelSize * pixelSize);
        
        // Create particle analyzer
        ResultsTable rt = new ResultsTable();
        ParticleAnalyzer pa = new ParticleAnalyzer(
            ParticleAnalyzer.ADD_TO_MANAGER,
            ParticleAnalyzer.AREA + ParticleAnalyzer.PERIMETER + ParticleAnalyzer.CIRCULARITY,
            rt,
            minAreaPixels,
            maxAreaPixels,
            minCircularity, maxCircularity
        );
        
        ImagePlus maskImp = new ImagePlus("mask", segmented);
        pa.analyze(maskImp);
        
        Roi[] detectedRois = roiManager.getRoisAsArray();
        
        Overlay overlay = new Overlay();
        
        // Display probe region (only red outline, no fill, not added to ROI Manager)
        if (userROI != null) {
            Roi probeRoiCopy = (Roi) userROI.clone();
            probeRoiCopy.setStrokeColor(Color.RED);
            probeRoiCopy.setStrokeWidth(2);
            probeRoiCopy.setFillColor(null);
            overlay.add(probeRoiCopy);
        }
        
        int validCellCount = 0;
        int particleCount = rt.size();
        
        // Clear ROI Manager, prepare to add only eligible cells
        roiManager.reset();
        
        for (int i = 0; i < detectedRois.length; i++) {
            if (i >= particleCount) continue;
            
            Roi roi = detectedRois[i];
            
            Rectangle bounds = roi.getBounds();
            double centerX = bounds.x + bounds.width / 2.0;
            double centerY = bounds.y + bounds.height / 2.0;
            Point2D.Double cellCentroid = new Point2D.Double(centerX, centerY);
            
            double deviation = calculateDeviation(roi, cellCentroid, probePolygon, 
                                                probeCentroid, probeToEdgeDistances);
            
            double areaPixels = rt.getValue("Area", i);
            double perimeterPixels = rt.getValue("Perim.", i);
            double circularity = rt.getValue("Circ.", i);
            
            double areaReal = areaPixels * pixelSize * pixelSize;
            double perimeterReal = perimeterPixels * pixelSize;
            
            if (deviation <= distanceDeviation) {
                // Cells near the pentagon - purple fill, add to ROI Manager
                validCellCount++;
                roi.setStrokeColor(Color.MAGENTA);
                roi.setStrokeWidth(2);
                roi.setFillColor(new Color(128, 0, 128, 100));
                
                String roiName = "1f-cell" + validCellCount;
                roi.setName(roiName);
                
                ParticleData particle = new ParticleData(validCellCount, 
                    areaReal, perimeterReal, circularity, deviation, cellCentroid, roi);
                particles.put(validCellCount, particle);
                
                // Add to ROI Manager
                roiManager.addRoi(roi);
                
                TextRoi label = new TextRoi(
                    bounds.x + 5,
                    bounds.y - 10,
                    roiName,
                    new Font("Arial", Font.BOLD, 14)
                );
                label.setStrokeColor(Color.MAGENTA);
                overlay.add(label);
                
                log(String.format("Cell %d (near): deviation=%.1f%%, area=%.1f μm²", 
                    validCellCount, deviation, areaReal));
            } else {
                // Cells far from pentagon - yellow fill, not added to ROI Manager
                roi.setStrokeColor(Color.YELLOW);
                roi.setStrokeWidth(1);
                roi.setFillColor(new Color(255, 255, 0, 50));
                overlay.add(roi);
            }
        }
        
        previewWindow.setOverlay(overlay);
        previewWindow.show();
        
        roiManager.runCommand("Show All");
        
        log(String.format("Total cells detected: %d, Eligible cells: %d", particleCount, validCellCount));
    }
    
    private void trackAllFrames() {
        int nFrames = currentImp.getStackSize();
        
        IJ.showProgress(0);
        IJ.showStatus("Tracking particles...");
        
        for (int frame = 2; frame <= nFrames; frame++) {
            IJ.showProgress(frame, nFrames);
            
            ImageProcessor ip = currentImp.getStack().getProcessor(frame).duplicate();
            ImageProcessor segmented = segmentCells(ip);
            
            for (Map.Entry<Integer, ParticleData> entry : particles.entrySet()) {
                int particleId = entry.getKey();
                ParticleData particle = entry.getValue();
                
                if (particle.dropFrame != -1) continue;
                
                Rectangle bounds = particle.initialRoi.getBounds();
                
                int expand = 5;
                Rectangle searchBounds = new Rectangle(
                    Math.max(0, bounds.x - expand),
                    Math.max(0, bounds.y - expand),
                    Math.min(ip.getWidth() - bounds.x + expand, bounds.width + 2 * expand),
                    Math.min(ip.getHeight() - bounds.y + expand, bounds.height + 2 * expand)
                );
                
                if (!isParticlePresent(segmented, searchBounds)) {
                    particle.dropFrame = frame;
                    particle.dropTime = (frame - 1) * frameInterval / 1000.0;
                    final int finalFrame = frame;
                    final int finalId = particle.id;
                    SwingUtilities.invokeLater(new Runnable() {
                        @Override
                        public void run() {
                            log("Cell " + finalId + " detached at frame " + finalFrame);
                        }
                    });
                }
            }
        }
        
        IJ.showProgress(1.0);
        IJ.showStatus("Tracking completed");
    }
    
    private boolean isParticlePresent(ImageProcessor segmented, Rectangle searchBounds) {
        int whitePixels = 0;
        int totalPixels = 0;
        
        int startX = Math.max(0, searchBounds.x);
        int startY = Math.max(0, searchBounds.y);
        int endX = Math.min(segmented.getWidth(), searchBounds.x + searchBounds.width);
        int endY = Math.min(segmented.getHeight(), searchBounds.y + searchBounds.height);
        
        for (int y = startY; y < endY; y++) {
            for (int x = startX; x < endX; x++) {
                totalPixels++;
                if (segmented.getPixel(x, y) == 255) {
                    whitePixels++;
                }
            }
        }
        
        if (totalPixels == 0) return false;
        
        return (whitePixels * 100.0 / totalPixels) > 20;
    }
    
    private void calculateMaxDropTime() {
        maxDropTime = 0;
        for (ParticleData particle : particles.values()) {
            if (particle.dropTime > maxDropTime) {
                maxDropTime = particle.dropTime;
            }
        }
        // If no cells detached, set a default value
        if (maxDropTime <= 0) {
            maxDropTime = 1.0;
        }
        log("Maximum detachment time: " + String.format("%.2f seconds", maxDropTime));
    }
    
    /* ==================== Pseudo-color function ==================== */
    private void showPseudoColorFrame() {
        if (particles.isEmpty()) {
            return;
        }
        
        // If window already displayed, bring to front
        if (pseudoFrame != null && pseudoFrame.isVisible()) {
            pseudoFrame.toFront();
            return;
        }
        
        // Close old pseudo-color window
        if (pseudoColorWindow != null) {
            pseudoColorWindow.close();
        }
        
        // Initialize current thresholds
        currentMinThreshold = 0;
        currentMaxThreshold = maxDropTime * detachmentConstant;
        currentLutType = "Jet";
        
        // Create independent JFrame window
        pseudoFrame = new JFrame("Pseudo-color Settings");
        pseudoFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        pseudoFrame.setLayout(new BorderLayout());
        pseudoFrame.setSize(450, 350);
        pseudoFrame.setLocationRelativeTo(mainFrame);
        
        // Add window close listener
        pseudoFrame.addWindowListener(new java.awt.event.WindowAdapter() {
            @Override
            public void windowClosing(java.awt.event.WindowEvent e) {
                pseudoFrame = null;
            }
        });
        
        JPanel mainPanel = new JPanel();
        mainPanel.setLayout(new BoxLayout(mainPanel, BoxLayout.Y_AXIS));
        mainPanel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
        
        // LUT Preview label
        JLabel lutLabel = new JLabel("LUT Preview");
        lutLabel.setAlignmentX(Component.CENTER_ALIGNMENT);
        lutLabel.setFont(new Font("SansSerif", Font.BOLD, 12));
        
        // Pseudo-color bar preview (horizontal) - will stretch with window
        lutPreviewPanel = createLUTPreviewPanel(currentLutType);
        lutPreviewPanel.setBorder(BorderFactory.createLineBorder(Color.BLACK));
        lutPreviewPanel.setPreferredSize(new Dimension(400, 60));
        lutPreviewPanel.setMinimumSize(new Dimension(200, 40));
        lutPreviewPanel.setMaximumSize(new Dimension(Short.MAX_VALUE, 80));
        
        // Pseudo-color type selection - 移除Fire和Ice
        JPanel lutPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        lutPanel.add(new JLabel("LUT Type:"));
        String[] lutTypes = {"Jet", "Spectrum", "Red", "Green", "Blue"};
        JComboBox<String> lutCombo = new JComboBox<>(lutTypes);
        lutCombo.setSelectedIndex(0);
        lutCombo.addActionListener(e -> {
            String selected = (String) lutCombo.getSelectedItem();
            currentLutType = selected;
            lutPreviewPanel.repaint();
        });
        lutPanel.add(lutCombo);
        
        // Detachment constant setting
        JPanel constantPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        constantPanel.add(new JLabel("Detachment Constant C:"));
        JTextField constantField = new JTextField(String.format("%.2f", detachmentConstant), 10);
        constantPanel.add(constantField);
        
        // Threshold range settings
        JPanel thresholdPanel = new JPanel(new GridLayout(2, 2, 5, 5));
        thresholdPanel.setBorder(BorderFactory.createTitledBorder("Threshold Range"));
        
        thresholdPanel.add(new JLabel("Min Threshold (s):"));
        JTextField minField = new JTextField(String.format("%.2f", currentMinThreshold), 10);
        thresholdPanel.add(minField);
        
        thresholdPanel.add(new JLabel("Max Threshold (s):"));
        JTextField maxField = new JTextField(String.format("%.2f", currentMaxThreshold), 10);
        thresholdPanel.add(maxField);
        
        // Button panel - only Confirm button
        JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.CENTER));
        JButton okButton = new JButton("Confirm");
        
        okButton.addActionListener(e -> {
            try {
                double minThreshold = Double.parseDouble(minField.getText());
                double maxThreshold = Double.parseDouble(maxField.getText());
                double constant = Double.parseDouble(constantField.getText());
                String lutType = (String) lutCombo.getSelectedItem();
                
                // Close old pseudo-color window
                if (pseudoColorWindow != null) {
                    pseudoColorWindow.close();
                }
                
                // Update detachment constant
                detachmentConstant = constant;
                
                // Generate new pseudo-color image with scaled times
                generatePseudoColorImage(lutType, minThreshold, maxThreshold);
                
                // Update current settings
                currentMinThreshold = minThreshold;
                currentMaxThreshold = maxThreshold;
                currentLutType = lutType;
                
                // Update preview
                lutPreviewPanel.repaint();
                
            } catch (NumberFormatException ex) {
                IJ.showMessage("Please enter valid numbers");
            }
        });
        
        buttonPanel.add(okButton);
        
        // Assemble panels
        mainPanel.add(lutLabel);
        mainPanel.add(Box.createVerticalStrut(5));
        mainPanel.add(lutPreviewPanel);
        mainPanel.add(Box.createVerticalStrut(10));
        mainPanel.add(lutPanel);
        mainPanel.add(Box.createVerticalStrut(10));
        mainPanel.add(constantPanel);
        mainPanel.add(Box.createVerticalStrut(10));
        mainPanel.add(thresholdPanel);
        mainPanel.add(Box.createVerticalStrut(10));
        mainPanel.add(buttonPanel);
        
        pseudoFrame.add(mainPanel, BorderLayout.CENTER);
        pseudoFrame.setVisible(true);
        
        // Generate initial pseudo-color image
        generatePseudoColorImage(currentLutType, currentMinThreshold, currentMaxThreshold);
    }
    
    /**
     * Create LUT preview panel
     */
    private JPanel createLUTPreviewPanel(String lutType) {
        JPanel panel = new JPanel() {
            @Override
            protected void paintComponent(Graphics g) {
                super.paintComponent(g);
                int width = getWidth();
                int height = getHeight();
                
                if (width <= 0 || height <= 0) return;
                
                // Get the current LUT type from the class variable
                String currentType = currentLutType;
                
                for (int x = 0; x < width; x++) {
                    int value = (int)((double)x / width * 255);
                    Color color = getColorFromLUT(currentType, value);
                    g.setColor(color);
                    g.drawLine(x, 0, x, height - 15);
                }
                
                // Draw border
                g.setColor(Color.BLACK);
                g.drawRect(0, 0, width - 1, height - 16);
                
                // Label thresholds at bottom
                g.setColor(Color.BLACK);
                g.setFont(new Font("SansSerif", Font.PLAIN, 10));
                g.drawString(String.format("%.1f", currentMinThreshold), 5, height - 2);
                g.drawString(String.format("%.1f", currentMaxThreshold), width - 50, height - 2);
            }
        };
        panel.setPreferredSize(new Dimension(400, 70));
        panel.setBackground(Color.WHITE);
        return panel;
    }
    
    private void generatePseudoColorImage(String lutType, double minThreshold, double maxThreshold) {
        int width = currentImp.getWidth();
        int height = currentImp.getHeight();
        
        // Create RGB image with white background
        ColorProcessor cp = new ColorProcessor(width, height);
        // Set white background (RGB: 255,255,255)
        int white = 0xFFFFFF;
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                cp.set(x, y, white);
            }
        }
        
        ImagePlus pseudoImp = new ImagePlus("Detachment Time Pseudo-color Map", cp);
        
        // Create overlay for contours
        Overlay overlay = new Overlay();
        
        // Get ROIs from ROI Manager
        RoiManager roiManager = RoiManager.getInstance();
        if (roiManager != null) {
            Roi[] rois = roiManager.getRoisAsArray();
            
            // Fill color for each ROI
            for (Roi roi : rois) {
                String roiName = roi.getName();
                if (roiName != null && roiName.startsWith("1f-cell")) {
                    // Parse cell ID
                    int cellId = Integer.parseInt(roiName.substring(7));
                    
                    // Find corresponding particle data
                    for (ParticleData particle : particles.values()) {
                        if (particle.id == cellId) {
                            // Calculate scaled time: dropTime * C
                            double dropTime = particle.dropTime;
                            if (dropTime == -1) {
                                dropTime = maxThreshold / detachmentConstant;
                            }
                            double scaledTime = dropTime * detachmentConstant;
                            
                            // Determine color value (0-255) based on scaled time
                            int colorValue;
                            if (maxThreshold <= minThreshold) {
                                colorValue = 255;
                            } else {
                                double ratio = (scaledTime - minThreshold) / (maxThreshold - minThreshold);
                                ratio = Math.max(0, Math.min(1, ratio));
                                colorValue = (int) (ratio * 255);
                            }
                            
                            // Create fill ROI
                            Roi fillRoi = (Roi) roi.clone();
                            
                            // Get color based on LUT type
                            Color fillColor = getColorFromLUT(lutType, colorValue);
                            
                            fillRoi.setFillColor(fillColor);
                            fillRoi.setStrokeColor(null);
                            overlay.add(fillRoi);
                            
                            break;
                        }
                    }
                }
            }
        }
        
        // Add probe outline (black)
        if (userROI != null) {
            Roi probeRoiCopy = (Roi) userROI.clone();
            probeRoiCopy.setStrokeColor(Color.BLACK);
            probeRoiCopy.setStrokeWidth(2);
            probeRoiCopy.setFillColor(null);
            overlay.add(probeRoiCopy);
        }
        
        pseudoImp.setOverlay(overlay);
        pseudoImp.show();
        
        pseudoColorWindow = pseudoImp;
    }
    
    /**
     * Get color based on LUT type and gray value
     */
    private Color getColorFromLUT(String lutType, int value) {
        float ratio = value / 255.0f;
        
        switch (lutType) {
            case "Jet":
                // Jet with proportions: red-orange-yellow (40%), green (20%), cyan-blue (40%)
                if (ratio < 0.4) {
                    // Red to orange to yellow (0-0.4)
                    if (ratio < 0.2) {
                        // Red to orange (0-0.2)
                        return new Color(255, (int)(ratio * 5 * 255), 0);
                    } else {
                        // Orange to yellow (0.2-0.4)
                        return new Color(255, 255, 0);
                    }
                } else if (ratio < 0.6) {
                    // Green (0.4-0.6)
                    float greenRatio = (ratio - 0.4f) * 5f;
                    return new Color((int)((1 - greenRatio) * 255), 255, 0);
                } else {
                    // Cyan to blue (0.6-1.0)
                    if (ratio < 0.8) {
                        // Cyan (0.6-0.8)
                        float cyanRatio = (ratio - 0.6f) * 5f;
                        return new Color(0, 255, (int)(cyanRatio * 255));
                    } else {
                        // Blue (0.8-1.0)
                        float blueRatio = (ratio - 0.8f) * 5f;
                        return new Color(0, (int)((1 - blueRatio) * 255), 255);
                    }
                }
                
            case "Red":
                return new Color(value, 0, 0);
                
            case "Green":
                return new Color(0, value, 0);
                
            case "Blue":
                return new Color(0, 0, value);
                
            case "Spectrum":
            default:
                // Spectrum: red->orange->yellow->green->cyan->blue->purple
                int r, g, b;
                if (ratio < 0.2) {
                    r = 255;
                    g = (int)(ratio * 5 * 255);
                    b = 0;
                } else if (ratio < 0.4) {
                    r = (int)((1 - (ratio - 0.2) * 5) * 255);
                    g = 255;
                    b = 0;
                } else if (ratio < 0.6) {
                    r = 0;
                    g = 255;
                    b = (int)((ratio - 0.4) * 5 * 255);
                } else if (ratio < 0.8) {
                    r = 0;
                    g = (int)((1 - (ratio - 0.6) * 5) * 255);
                    b = 255;
                } else {
                    r = (int)((ratio - 0.8) * 5 * 255);
                    g = 0;
                    b = 255;
                }
                return new Color(r, g, b);
        }
    }
    
    /* ==================== 新策略：Results display method ==================== */
    private void displayResults() {
        // 新策略：每次使用不同的标题，避免窗口复用问题
        String uniqueTitle = baseResultsTitle + " - Analysis #" + analysisCounter;
        
        // 创建新的结果表格
        ResultsTable rt = new ResultsTable();
        
        for (ParticleData particle : particles.values()) {
            rt.incrementCounter();
            rt.addValue("Cell ID", "cell" + particle.id);
            rt.addValue("Initial Area (μm²)", String.format("%.2f", particle.initialArea));
            rt.addValue("Initial Perimeter (μm)", String.format("%.2f", particle.initialPerimeter));
            rt.addValue("Initial Circularity", String.format("%.3f", particle.initialCircularity));
            rt.addValue("Distance Deviation (%)", String.format("%.1f", particle.distanceDeviation));
            
            if (particle.dropTime == -1) {
                rt.addValue("Detachment Time (s)", "Not detached");
                rt.addValue("Detachment Frame", "Not detached");
            } else {
                rt.addValue("Detachment Time (s)", String.format("%.2f", particle.dropTime));
                rt.addValue("Detachment Frame", particle.dropFrame);
            }
        }
        
        // 显示新表格
        rt.show(uniqueTitle);
        
        // 不保存引用，让ImageJ自己管理
    }
    
    private void showSummary() {
        int totalCells = particles.size();
        int droppedCells = 0;
        double totalDropTime = 0;
        
        for (ParticleData particle : particles.values()) {
            if (particle.dropTime != -1) {
                droppedCells++;
                totalDropTime += particle.dropTime;
            }
        }
        
        double dropRatio = totalCells > 0 ? (double)droppedCells / totalCells * 100 : 0;
        double avgDropTime = droppedCells > 0 ? totalDropTime / droppedCells : 0;
        
        // Close previous summary dialog
        if (summaryDialog != null) {
            summaryDialog.dispose();
        }
        
        summaryDialog = new JDialog(mainFrame, "Analysis Summary", false);
        summaryDialog.setLayout(new BorderLayout());
        summaryDialog.setSize(450, 350);
        summaryDialog.setLocationRelativeTo(mainFrame);
        summaryDialog.setResizable(true);
        
        JTextArea textArea = new JTextArea();
        textArea.setEditable(false);
        textArea.setFont(new Font("SansSerif", Font.PLAIN, 14));
        
        StringBuilder summary = new StringBuilder();
        summary.append("\n");
        summary.append("╔════════════════════════════════════════╗\n");
        summary.append("║    Probe Region Cell Detachment Summary    ║\n");
        summary.append("╚════════════════════════════════════════╝\n\n");
        
        summary.append("【Analysis Parameters】\n");
        summary.append(String.format("  Scale: 1 pixel = %.4f μm\n", pixelSize));
        summary.append(String.format("  Probe region: Set (pentagon)\n"));
        summary.append(String.format("  Cell distance deviation threshold: %.1f%%\n", distanceDeviation));
        summary.append(String.format("  Area range: %.1f - %.1f μm²\n", minArea, maxArea));
        summary.append(String.format("  Circularity range: %.2f - %.2f\n", minCircularity, maxCircularity));
        summary.append(String.format("  Frame interval: %.1f ms\n\n", frameInterval));
        
        summary.append("【Statistics】\n");
        summary.append(String.format("  Total cells: %d\n", totalCells));
        summary.append(String.format("  Detached cells: %d\n", droppedCells));
        summary.append(String.format("  Detachment ratio: %.1f%%\n", dropRatio));
        summary.append(String.format("  Average detachment time: %.2f s\n", avgDropTime));
        summary.append(String.format("  Maximum detachment time: %.2f s\n\n", maxDropTime));
        
        summary.append("Analysis complete! See table for full results.\n");
        
        textArea.setText(summary.toString());
        textArea.setMargin(new Insets(10, 10, 10, 10));
        
        JScrollPane scrollPane = new JScrollPane(textArea);
        summaryDialog.add(scrollPane, BorderLayout.CENTER);
        
        JButton closeButton = new JButton("Close");
        closeButton.addActionListener(e -> summaryDialog.dispose());
        JPanel buttonPanel = new JPanel();
        buttonPanel.add(closeButton);
        summaryDialog.add(buttonPanel, BorderLayout.SOUTH);
        
        summaryDialog.setVisible(true);
    }
}
