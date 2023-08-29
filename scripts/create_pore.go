package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"strconv"
	"strings"
)

func boxVectors(inputFilePath string) ([3]float64, error) {
	// Open the file
	file, err := os.Open(inputFilePath)
	if err != nil {
		return [3]float64{}, err
	}
	defer file.Close()

	// Create a scanner to read the file line by line
	scanner := bufio.NewScanner(file)

	var lastLine string

	// Read the file line by line
	for scanner.Scan() {
		lastLine = scanner.Text()
	}

	if err := scanner.Err(); err != nil {
		return [3]float64{}, err
	}

	// Split the last line into float values
	floatStrings := strings.Fields(lastLine)
	if len(floatStrings) != 3 {
		return [3]float64{}, fmt.Errorf("The last line does not contain three float values")
	}

	// Convert float string values to actual float variables
	var boxDims [3]float64
	for i, floatStr := range floatStrings {
		float, err := strconv.ParseFloat(floatStr, 64)
		if err != nil {
			return [3]float64{}, err
		}
		boxDims[i] = float
	}

	return boxDims, nil
}

func shouldRemove(moleculeName string, molecules []string) bool {
	for _, molecule := range molecules {
		if strings.HasSuffix(moleculeName, molecule) {
			return true
		}
	}
	return false
}

func stringInSlice(a string, list []string) bool {
    for _, b := range list {
        if b == a {
            return true
        }
    }
    return false
}


func parseGroLine(line string) (atomName, atomType string, x, y, z float64) {
	atomName = line[0:5]
	atomType = line[5:10]
	x, _ = strconv.ParseFloat(strings.TrimSpace(line[20:28]), 64)
	y, _ = strconv.ParseFloat(strings.TrimSpace(line[28:36]), 64)
	z, _ = strconv.ParseFloat(strings.TrimSpace(line[36:44]), 64)
	return atomName, atomType,  x, y, z
}

func formPore(inputFilePath string, outputFilePath string, axis string, center1 float64, center2 float64, radius float64, molecules []string) (int, error) {

	// Read the input GROMACS gro file
	file, err := os.Open(inputFilePath)
	if err != nil {
		return 0, err
	}
	defer file.Close()

	var lines []string
	var removedMolecules int
	var molsToRemove []string
	var currentMoleculeLines []string
	var firstLine, lastLine string

	isWithinCircle := func(dim1, dim2 float64, boxDims [3]float64) bool {
        dx := dim1 - center1
        dy := dim2 - center2

        // Handle periodic boundary conditions
        if dx > boxDims[0]/2 {
            dx -= boxDims[0]
        } else if dx < -boxDims[0]/2 {
            dx += boxDims[0]
        }

        if dy > boxDims[1]/2 {
            dy -= boxDims[1]
        } else if dy < -boxDims[1]/2 {
            dy += boxDims[1]
        }

        return dx*dx + dy*dy <= radius*radius
    }

	boxDims, _ := boxVectors(inputFilePath)

	fmt.Printf("Drilling a pore in the %v axis at %.2f, %.2f, with radius %.1f\n", axis, center1, center2, radius)

	scanner := bufio.NewScanner(file)
	lineCount := 0
	for scanner.Scan() {
		line := scanner.Text()

		// Store the first and last lines
		if lineCount == 0 {
			firstLine = line
		}
		lastLine = line

		if lineCount <= 1 || (len(line) != 44 && len(line) != 68){  // skip first two lines or any line that doesn't fit gromacs file req (lazy, but allows skipping last line or writing comments idk)
			lineCount++
			continue
		} else  {
			atomName, atomType, x, y, z := parseGroLine(line)
			var resid string
			resid = atomName + atomType

			var dim1, dim2 float64

			if axis == "x" {
				dim1  = y
				dim2  = z
			} else if axis == "y" {
				dim1  = x
				dim2  = z
			} else if axis == "z" {
				dim1  = x
				dim2  = y
			}

			if shouldRemove(resid, molecules) && isWithinCircle(dim1, dim2, boxDims) {
				if stringInSlice(resid, molsToRemove) == false {
					molsToRemove = append(molsToRemove, resid)
					fmt.Printf("One or more atoms of molecule %s falls within pore radius, molecule will be deleted.\n", resid)
					// fmt.Printf("Removing molecule %s\n", resid)
					removedMolecules++
				} 
			}

			currentMoleculeLines = append(currentMoleculeLines, line)

			lineCount++

		}
	}

	if err := scanner.Err(); err != nil {
		return 0, err
	}

    // Prepare the lines to write back to the output file
    lines = append(lines, firstLine)                     // Add the first line from the original input file
    lines = append(lines, "")                            // Placeholder for the line count
	for _, line := range currentMoleculeLines {   		 // Add the current molecule lines
		if !stringInSlice((line[:10]), molsToRemove) {  
			lines = append(lines, line)
		}
	}
    // lines = append(lines, currentMoleculeLines...)       // Add the current molecule lines
    lines = append(lines, lastLine)                      // Add the last line from the original input file

    // Calculate the line count after processing
    lines[1] = fmt.Sprintf("%d", len(lines)-3)

    // Write the updated gro file
    // outputFilePath := strings.TrimSuffix(inputFilePath, ".gro") + "_modified.gro"
    outputFile, err := os.Create(outputFilePath)
    if err != nil {
        return 0, err
    }
    defer outputFile.Close()


	for _, line := range lines {
		fmt.Fprintln(outputFile, line)
	}
	

    return removedMolecules, nil
}

func contains(slice []string, element string) bool {
	for _, item := range slice {
		if item == element {
			return true
		}
	}
	return false
}


func main() {
	inputFilePath := flag.String("c", "", "Path to the input file")
	outputFilePath := flag.String("o", strings.TrimSuffix(*inputFilePath, ".gro") + "_modified.gro", "Path to the output file")
	axis := flag.String("axis", "z", "Axis to print values for (x, y, or z)")
	poreRadius := flag.Float64("r", 1, "Pore radius as an integer")
	poreCenter := flag.String("poreCenter", "", "Pore center coordinates (two floats separated by a comma)")
	flag.Parse()

	if *inputFilePath == "" {
		fmt.Println("Usage: program_name -c input_file_path -axis x|y|z -r pore_radius [-poreCenter center1,center2]")
		return
	}


	var center1, center2 float64

	if *poreCenter != "" {
		// Parse and split the comma-separated pore center values
		centerValues := strings.Split(*poreCenter, ",")
		if len(centerValues) != 2 {
			fmt.Println("Invalid pore center coordinates.")
			return
		}

		// Convert string values to float64
		var err error
		center1, err = strconv.ParseFloat(centerValues[0], 64)
		if err != nil {
			fmt.Println("Error parsing center1:", err)
			return
		}
		center2, err = strconv.ParseFloat(centerValues[1], 64)
		if err != nil {
			fmt.Println("Error parsing center2:", err)
			return
		}
	} else {
		boxDims, err := boxVectors(*inputFilePath)
		if err != nil {
			fmt.Println("Error:", err)
			return
		}

		switch *axis {
		case "x":
			center1 = boxDims[1]/2
			center2 = boxDims[2]/2
		case "y":
			center1 = boxDims[0]/2
			center2 = boxDims[2]/2
		case "z":
			center1 = boxDims[0]/2
			center2 = boxDims[1]/2
		default:
			fmt.Println("Invalid axis value. Choose x, y, or z.")
			return
		}
	}

	removedMolecules, err := formPore(*inputFilePath, *outputFilePath, *axis, center1, center2, *poreRadius, []string{""})
	if err != nil {
		fmt.Println("Error:", err)
	} else {
		fmt.Printf("Removed %d molecules\n", removedMolecules)
	}

}