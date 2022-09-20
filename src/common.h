#pragma once

inline double piano_key_frequency(uint8_t note_index)
{
    return pow(pow(2.0, 1.0 / 12.0), note_index);
}

inline uint8_t keyboard_to_piano_key(char key)
{
    static const uint8_t keyboard_piano_key_mapping[] = {
        7, 4, 2, 9, 18,
        10, 11, 12, 23, 13,
        14, 15, 6, 5, 24,
        25, 16, 19, 8, 20,
        22, 3, 17, 1, 21, 0
    };

    return keyboard_piano_key_mapping[key - 'a'];
}

inline double keyboard_to_piano_key_frequency(char key)
{
    return piano_key_frequency(keyboard_to_piano_key(key));
}
