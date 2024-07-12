// grab cookie
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
// fetchapi
async function fetchAPI(url, formData, csrftoken) {
  const response = await fetch(url, {
    method: "POST",
    body: formData,
    headers: {
      "X-CSRFToken": csrftoken,
    },
  });

  if (!response.ok) {
    toggleLoading(false);
    throw new Error("Network response was not ok");
  }
  return response.json();
}

// fetching image
function loadImage(folder, filename, containerId) {
  console.log(containerId);
  const csrftoken = getCookie("csrftoken");
  fetch("/get_image/", {
    method: "POST",
    headers: {
      "Content-Type": "application/x-www-form-urlencoded",
      "X-CSRFToken": csrftoken,
    },
    body: `folder=${encodeURIComponent(folder)}&filename=${encodeURIComponent(
      filename
    )}`,
  })
    .then((response) => response.json())
    .then((data) => {
      let container = document.getElementById(containerId);
      if (container) {
        if (data.image_path) {
          if (document.getElementById(filename)) {
            document.getElementById(filename).remove();
          }
          let img = document.createElement("img");
          img.src = data.image_path;
          img.alt = filename;
          img.id = filename;
          img.classList.add("h-auto", "max-w-full", "mx-auto");
          container.appendChild(img);
        } else {
          container.innerHTML = "No image found. Might be something got wrong.";
        }
      } else {
        console.error("Container not found:", containerId);
      }
    })
    .catch((error) => {
      console.error("Error:", error);
      let container = document.getElementById(containerId);
      if (container) {
        container.innerHTML = "No image found. Might be something got wrong.";
      }
    });
}

// // 使用示例
// document.addEventListener("DOMContentLoaded", function () {
//   loadImage('clustering_result', 'heatmap.png', 'heatmap-container');

//   // 你可以多次調用 loadImage 來加載不同的圖片
//   // loadImage('another_folder', 'another_image.jpg');
// });
