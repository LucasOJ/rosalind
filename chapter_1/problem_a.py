# https://rosalind.info/problems/ba1a/

def pattern_count(text:str, pattern:str) -> int:
    text_len = len(text)
    pattern_len = len(pattern)

    count = 0

    for i in range(0, text_len - pattern_len):
        if text[i:i+pattern_len] == pattern:
            count += 1

    return count
