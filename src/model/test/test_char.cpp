#include <iostream>
#include <cstring>

const char** stringToConstCharPtrArray(const std::string& str) {
    // Allocate memory for a const char* array with one element
    const char** strArray = new const char*[1];

    // Create a copy of the string as a const char* and store it in the array
    strArray[0] = strdup(str.c_str());

    return strArray;
}

int main() {
    std::string myString = "Hello, World!";
    
    const char** strArray = stringToConstCharPtrArray(myString);

    // Print the value in the const char* array
    std::cout << "Converted String: " << strArray[0] << std::endl;

    // Don't forget to free memory when done to avoid memory leaks
    free((void*)strArray[0]); // Free the duplicated string
    delete[] strArray; // Free the array itself

    return 0;
}
