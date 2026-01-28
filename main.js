function copyToClipboard(text) {
    navigator.clipboard.writeText(text).then(() => {
        console.log('Text copied to clipboard');
        // Optional: Add a visual indicator to show success
        const buttons = document.querySelectorAll('.copy-btn');
        buttons.forEach(button => {
            const originalIcon = button.innerHTML;
            if (button.nextElementSibling.innerText.includes(text) || button.parentElement.querySelector('pre code').innerText.includes(text)) {
                button.innerHTML = '<i class="fa-solid fa-check"></i> Copied!';
                setTimeout(() => {
                    button.innerHTML = originalIcon;
                }, 2000);
            }
        });
    }).catch(err => {
        console.error('Failed to copy text: ', err);
    });
}

// Smooth scrolling for navigation links
document.querySelectorAll('a[href^="#"]').forEach(anchor => {
    anchor.addEventListener('click', function (e) {
        e.preventDefault();

        document.querySelector(this.getAttribute('href')).scrollIntoView({
            behavior: 'smooth'
        });
    });
});