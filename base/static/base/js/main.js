function getCookie(name) {
  let cookieValue = null;
  if (document.cookie && document.cookie !== "") {
    const cookies = document.cookie.split(";");
    for (let i = 0; i < cookies.length; i++) {
      const cookie = cookies[i].trim();
      if (cookie.substring(0, name.length + 1) === name + "=") {
        cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
        break;
      }
    }
  }
  return cookieValue;
}

let pathSegments = window.location.pathname.split("/");
let secondSegment = pathSegments[2] || "";

let navContents = document.getElementsByName("nav-content");

for (let i = 0; i < navContents.length; i++) {
  let element = navContents[i];
  let content = element.textContent.replace(/\s+/g, "").toLowerCase();

  if (content === secondSegment) {
    element.classList.add("text-sky-700");
    break;
  }
  element.classList.add("text-sky-700");
}
