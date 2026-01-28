// Copy to Clipboard Functionality
function copyToClipboard(text) {
    navigator.clipboard.writeText(text).then(() => {
        // Find the button that was clicked (approximate logic since we pass text)
        // In a real app, we would pass 'this' or the event.
        // For simplicity, we just show a generic toast or alert, or update all buttons temporarily.
        
        // Let's iterate all buttons and find which one might have triggered this or just show a visual feedback
        const buttons = document.querySelectorAll('.copy-btn');
        buttons.forEach(btn => {
            const originalHTML = btn.innerHTML;
            // Simplified: we just briefly change icon
            if(btn.onclick.toString().includes(text.substring(0, 10))) {
                btn.innerHTML = '<i class="fa-solid fa-check"></i>';
                setTimeout(() => {
                    btn.innerHTML = '<i class="fa-regular fa-copy"></i>';
                }, 2000);
            }
        });
    }).catch(err => {
        console.error('Failed to copy: ', err);
    });
}

// Smooth Scroll for Anchors
document.querySelectorAll('a[href^="#"]').forEach(anchor => {
    anchor.addEventListener('click', function (e) {
        e.preventDefault();
        document.querySelector(this.getAttribute('href')).scrollIntoView({
            behavior: 'smooth'
        });
    });
});

// Simple Intersection Observer for Fade-in animations
const observerOptions = {
    threshold: 0.1
};

const observer = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
        if (entry.isIntersecting) {
            entry.target.style.opacity = '1';
            entry.target.style.transform = 'translateY(0)';
        }
    });
}, observerOptions);

document.querySelectorAll('.feature-card, .step').forEach(el => {
    el.style.opacity = '0';
    el.style.transform = 'translateY(20px)';
    el.style.transition = 'opacity 0.6s ease-out, transform 0.6s ease-out';
    observer.observe(el);
});
